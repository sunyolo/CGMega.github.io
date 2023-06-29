# Tutorial

This tutorial demonstrates how to use CGMega functions with a demo dataset. Once you are familiar with CGMega’s workflow, please replace the demo data with your own data to begin your analysis.

#What the tutorial covers

## 1. How to prepare input data

We recommend getting started with CGMega using the provided demo dataset. When you want to apply CGMega to your own multi-omics dataset, please refer to the following tutorials to learn how to prepare input data. Overall, the input data consists of two parts: the graph, constructed from PPI and the node feature including compressed Hi-C features, SNVs, CNVs frequencies and epigenetic densities.

#### PPI graph construction
First, we read the PPI data (from CPDB as default) and transform it into a graph through the following commands (in [data_preprocess_cv.py](https://github.com/NBStarry/CGMega/blob/main/data_preprocess_cv.py))
 ```
 # Load PPI data
 ppi_mat = get_ppi_mat(ppi, drop_rate=ppi_drop_rate, from_list=False, random_seed=random_seed, pan=pan) if ppi else None
 ```

 ```
 # Corresponding functions to read PPI data in different forms.
 def get_ppi_mat(ppi_name='CPDB', drop_rate=0.0, from_list=False, random_seed=42, reduce=False, pan=False):
    """
    Read PPI data from a csv file and construct PPI adjacent matrix.

    Parameters:
    ----------
    ppi_name:   str, {'CPDB', 'IRef', 'Multinet', 'PCNet', 'STRING'}, default='CPDB'. 
                Name of PPI network dataset. Corresponding matrix should be put into certain directory first.
    drop_rate:  float, default=0.0. 
                Drop rate for robustness study.
    from_list:  bool, default=False.
                Whether the PPI data is loaded from a preprocessed adjacency list instead of a matrix.
    random_seed:int, default=42.
    reduce:     bool, default=False.
                Whether to load ppi embedding got by Node2Vec.
    pan:        bool, default=False.
                Whether to use pan data in EMOGI.

    Returns:
    numpy.ndarray
        A ndarray(num_nodes, num_nodes) contains PPI adjacent matrix.
    """
    prefix = "data" if not pan else "pan_data"
    if reduce:
        ppi_dir = prefix + f"/{ppi_name}/N2V_ppi_embedding_15.csv"
        print(f"Loading PPI feature from {ppi_dir} ......")
        return read_table_to_np(ppi_dir)
    # Load PPI data from an edge list
    if from_list:
        ppi_dir = prefix + f"/{ppi_name}/{ppi_name}_edgelist.csv"
        print(f"Loading PPI matrix from {ppi_dir} ......")
        data = pd.read_csv(
            prefix + f"/{ppi_name}/{ppi_name}_edgelist.csv", sep='\t')
        # Load the gene names
        gene_list, gene_set = get_all_nodes(pan=pan)

        # Extract the edges that are also in the list of gene names
        if not pan:
            adj = [(row[1], row[2], row[3]) for row in data.itertuples()
                if row[1] in gene_set and row[2] in gene_set]
            conf = [row[4] for row in data.itertuples() if row[1]
                    in gene_set and row[2] in gene_set]
            if drop_rate:
                # Drop samples with stratification by confidence score
                adj, drop_adj = train_test_split(
                    adj, test_size=drop_rate, random_state=random_seed, stratify=conf)
            # Construct the adjacency matrix from the edges
            adj_matrix = pd.DataFrame(0, index=gene_list, columns=gene_list)
            for line in adj:
                adj_matrix.loc[line[0], line[1]] = line[2]
                adj_matrix.loc[line[1], line[0]] = line[2]
        else:
            adj = [(row[1], row[2]) for row in data.itertuples()
                if row[1] in gene_set and row[2] in gene_set]
            adj_matrix = pd.DataFrame(0, index=gene_list, columns=gene_list)
            for line in adj:
                adj_matrix.loc[line[0], line[1]] = 1
                adj_matrix.loc[line[1], line[0]] = 1
        adj_matrix.to_csv(prefix + f'/{ppi_name}/{ppi_name}_matrix.csv', sep='\t')
        data = adj_matrix.to_numpy()

        return data

    # Load PPI data from a matrix
    ppi_dir = prefix + f"/{ppi_name}/{ppi_name}_matrix.csv"
    print(f"Loading PPI matrix from {ppi_dir} ......")
    data = pd.read_csv(ppi_dir, sep="\t").to_numpy()[:, 1:]

    return data
 ```
 
 ```
 # Corresponding function to transform PPI data into a graph.
 def construct_edge(mat):
    """
    Construct edges from adjacent matrix.

    Parameters:
    ----------
    mat:    ndarray(num_nodes, num_nodes).
                PPI matrix from get_ppi_mat().

    Returns:
    edges:      list(num_edges, 2). 
    edge_dim:   int.
                Dim of edge features.
    val:        list(num_edges, ).
                Edge features(=[1] * num_edges in current version).
    """
    num_nodes = mat.shape[0]
    edges = []
    val = []
    for i in range(num_nodes):
        for j in range(num_nodes):
            if mat[i, j] > 0:
                edges.append([i, j])
                val.append(mat[i, j])

    edge_dim = 1
    edges = np.transpose(edges)
    val = np.reshape(val, (-1, edge_dim))

    return edges, edge_dim, val
 ```
 
 ```
 # use pyg
 def build_pyg_data(node_mat, node_lab, mat, pos):
    x = t.tensor(node_mat, dtype=torch.float)
    y = t.tensor(node_lab, dtype=torch.long)
    pos = t.tensor(pos, dtype=torch.int)
    edge_index, edge_dim, edge_feat = construct_edge(mat)
    edge_index = t.tensor(edge_index, dtype=torch.long)
    edge_feat = t.tensor(edge_feat, dtype=torch.float)
    data = Data(x=x.clone(), y=y.clone(), edge_index=edge_index,
                edge_attr=edge_feat, pos=pos, edge_dim=edge_dim)
    print(
        f"Number of edges: {data.num_edges}, Dimensionality of edge: {edge_dim},\nNubmer of nodes: {data.num_nodes}")

    return data
 ```
#### Hi-C data embedding

#### Other omics data

## 2. Model training

This section demonstrates how to train the GAT-based cancer gene prediction model.

- To reduce the number of parameters and make training feasible within time and resource constraints, the input graphs were sampled using neighbor sampler. The subgraphs included all first and second order neighbors for each node and training was performed on these subgraphs.
- The learning rate is increased linearly from 0 to 0.005 for the first 20% of the total iterations.
- warm-up strategy for learning rate is employed during the initial training phase.
- To prevent overfitting and over smoothing, an early stop strategy is adopted. If the model's performance on the validation set dose not improve for a consecutive 100 epochs, the training process stops.
- Dropout is used and the dropout rate is set to 0.1 for the attention mechanism and 0.4 for the other modules.
- Max pooling step size is 2. After pooling, the representation had 32 dimensions.

## 3. Interpretation

After prediction, you can do analyses as following to interpret your results.

#### 1. identify the gene module

#### 2. calculate the importance score

#### 3. calculate the Representative Feature

#### 4. explore the relationships between different gene modules

source: `{{ page.path }}`
