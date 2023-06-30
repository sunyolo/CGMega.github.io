# Tutorial

This tutorial demonstrates how to use CGMega functions with a demo dataset. Once you are familiar with CGMega’s workflow, please replace the demo data with your own data to begin your analysis.

<style>
pre {
  overflow-y: auto;
  max-height: 400px;
}
</style>

## 1. How to prepare input data

We recommend getting started with CGMega using the provided demo dataset. When you want to apply CGMega to your own multi-omics dataset, please refer to the following tutorials to learn how to prepare input data. Overall, the input data consists of two parts: the graph, constructed from PPI and the node feature including compressed Hi-C features, SNVs, CNVs frequencies and epigenetic densities.
 
 If you are unfamiliar with CGMega, you may start with our data used in the paper to save your time. For MCF7 cell line, K562 cell line and AML patients, the input data as well as their label information are uploaded [here](https://github.com/NBStarry/CGMega/tree/main/data). If you start with any one from these data, you can skip the step 1 about _How to prepare input data_.
 
 ```note
The labels should be collected yourself if you choose analyze your own data.
 ```
 
#### 1.1 Hi-C data embedding

We use SVD to condense the chromatin information in Hi-C data.

 ```note
Before SVD, the Hi-C data should go through: 

first, processing by [NeoLoopFinder](https://github.com/XiaoTaoWang/NeoLoopFinder) to remove the potential effects of structural variation; 

second, normalization using [ICE]( https://bitbucket.org/mirnylab/hiclib) correction to improve the data quality. If you are new to these two tools, please go through these document in advance.

[tutorial for NeoLoopFinder need link](./)

[tutorial for Hi-C normalization need link](./)
 ```

 ```css
 def get_hic_mat(data_dir='data/Breast_Cancer_Matrix', drop_rate=0.0, reduce='svd', reduce_dim=5, resolution='10KB', type='ice', norm='log'):
    """
    Read Hi-C matrix from a csv file and do dimensionality reduction or normalization. Corresponding matrix should be put into certain directory first.

    Parameters:
    ----------
    data_dir:   str, default='data/Breast_Cancer_Matrix'
    drop_rate:  float, default=0.0. 
                The proportion of entries to drop randomly for robustness study, set from 0.0 to 0.9. 
    reduce:     str, {'svd', 'svdr', 'nmf', 't-sne', 'umap', 'isomap', 'lle', False}, default='svd'. 
                Method for dimensionality reduction. Use False for no reduction.
    reduce_dim: int, default=5. 
                Dimensionality after reduction.
    resolution: str, {'5KB', '10KB', '25KB'}, default='10KB'.
                The resolution of Hi-C data.
    type:       str, {'ice', 'int'}, default='ice'.
                The type of Hi-C data.
    norm:       {'log', 'square', 'binary'}
                Whether to do min-max normalization for reduced Hi-C data.

    Returns:
    numpy.ndarray
        A matrix of shape (num_nodes, reduce_dim) if `reduce` is not False, or (num_nodes, num_nodes) otherwise.
    """
    DEFAULT_RESO = '10KB'
    cell_line = get_cell_line(data_dir)
    sample_rate = str(round(10 * (1-drop_rate)))

    def get_hic_dir():
        if type == 'ice':
            if sample_rate == '10' and resolution == DEFAULT_RESO:
                hic_dir = data_dir + cell_line + "_Adjacent_Matrix_Ice"
            else:
                hic_dir = data_dir + "/drop_hic_ice/" + resolution + '/' + \
                    sample_rate + '_' + cell_line[1:] + "_ICE_downsample.csv"
        elif type == 'int':
            hic_dir = data_dir + cell_line + "_Adjacent_Matrix"
        print(f"Loading Hi-C matrix from {hic_dir} ......")

        return hic_dir

    def normalize_hic_data(data, method, reduce):
        """
        Normalize the input Hi-C matrix.

        Parameters:
        ----------
        data: numpy.ndarray
              The input Hi-C matrix to be normalized.
        method:  str, {'log', 'square', 'binary'}
              The normalization method to use.

        Returns:
        -------
        numpy.ndarray
            The normalized Hi-C matrix.
        """
        UP_THRESHOLD = 10000
        DOWN_THRESHOLD = 5
        LOG_THRESHOLD = 0.7
        if method == 'log':
            data = np.log10(data + EPS)
            if reduce != False: return data
            else: 
                data[data<LOG_THRESHOLD] = 0
        elif method == 'square':
            # Clip the data to the range of [1, 10000] and apply a power-law transformation
            data[data > UP_THRESHOLD] = 10000
            data[data <= DOWN_THRESHOLD] = 0
            data = data ** 0.1
        elif method == 'binary':
            data = np.log10(data + EPS)
            data[data < LOG_THRESHOLD] = 0
            data[data > LOG_THRESHOLD] = 1
        else:
            raise ValueError(f"Invalid use value: {method}")

        return data

    if reduce == 'n2v':
        # Read pre-trained N2V embedding from a file
        hic_data = read_table_to_np(f'/N2V_embedding_{reduce_dim}.csv')
        return minmax(hic_data)

    reducer_dict = {
        'svd': TruncatedSVD(n_components=reduce_dim, algorithm="arpack"),
        'svdr': TruncatedSVD(n_components=reduce_dim),
        'nmf': NMF(n_components=reduce_dim, init='nndsvd', solver='mu',
                   beta_loss='frobenius', max_iter=10000, alpha=0.1, tol=1e-6, l1_ratio=1),
        't-sne': TSNE(n_components=reduce_dim,
                      learning_rate='auto', init='pca'),
        'umap': UMAP(n_components=reduce_dim),
        'isomap': Isomap(n_components=reduce_dim),
        'lle': LocallyLinearEmbedding(n_components=reduce_dim),
    }

    hic_data = read_table_to_np(
        get_hic_dir(), sep='\t', dtype=float, start_col=1)
    hic_data = normalize_hic_data(hic_data, method=norm, reduce=reduce)
    hic_data += 8 if reduce == 'nmf' else 0 

    reducer = reducer_dict[reduce] if reduce else None
    hic_data = reducer.fit_transform(hic_data) if reduce else hic_data

    return minmax(hic_data, axis= 1 if reduce else -1) 
 ```

#### 1.2 Other omics data

 ```
 code for other omics-data processing 
 ```

#### 1.3 PPI graph construction

Then,we read the PPI data (from CPDB as default) and transform it into a graph through the following commands (in [data_preprocess_cv.py](https://github.com/NBStarry/CGMega/blob/main/data_preprocess_cv.py)).

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

Now the input data is prepared. Let'a go for the model training!
 
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
