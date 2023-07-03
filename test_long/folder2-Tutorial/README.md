# Tutorial

This tutorial demonstrates how to use CGMega functions with a demo dataset (for MCF7 cell line).
Once you are familiar with CGMega’s workflow, please replace the demo data with your own data to begin your analysis.

<style>
pre {
  overflow-x: auto;
  overflow-y: auto;
  max-height: 300px;
}
</style>

## 1. How to prepare input data

We recommend getting started with CGMega using the provided demo dataset. When you want to apply CGMega to your own multi-omics dataset, please refer to the following tutorials to learn how to prepare input data.

Overall, the input data consists of two parts: the graph, constructed from PPI and the node feature including condensed Hi-C features, SNVs, CNVs frequencies and epigenetic densities.
 
 If you are unfamiliar with CGMega, you may start with our data used in the paper to save your time. For MCF7 cell line, K562 cell line and AML patients, the input data as well as their label information are uploaded [here](https://github.com/NBStarry/CGMega/tree/main/data). If you start with any one from these data, you can skip the _step 1_ about _How to prepare input data_.
 The following steps from 1.1~1.3 can be found in our source code [data_preprocess_cv.py](https://github.com/NBStarry/CGMega/blob/main/data_preprocess_cv.py)).
 ```note
The labels should be collected yourself if you choose analyze your own data.
 ```
 
 
#### 1.1 Hi-C data embedding

 ```note
Before SVD, the Hi-C data should go through: 
1. processing by [NeoLoopFinder](https://github.com/XiaoTaoWang/NeoLoopFinder) to remove the potential effects of structural variation; 
2. normalization using [ICE]( https://bitbucket.org/mirnylab/hiclib) correction to improve the data quality. 

If you are new to these two tools, please go through these document in advance.
[tutorial for NeoLoopFinder (need link)](./)
[tutorial for Hi-C normalization (need link)](./)
 ```
The parameters and main functions used in NeoLoopFinder are listed as below:

---
Parameters:
- input file format: .cool or .mcool
- resolution: 10Kb
- output file format: .matrix.
---
Functions:(availble in [batch_neoloop.sh]())
```
- calculate-cnv -H $hic.cool -g hg38 -e MboI --output ${hic}_10kb.cnv
- segment-cnv  --cnv-file ${hic}-10kb.cnv --binsize 10000 --output  ${hic}-10k.seg  --nproc 10  --ploidy 2
- cooler balance $hic.cool
- correct-cnv -H $hic.cool --cnv-file ${hic}-10k.seg  --nproc 10
- assemble-complexSVs -O ${hic}_10kb  -B $hic.sv  -H $hic.cool
- neoloop-caller -O $hic.neo-loops.txt  -H $hic.cool  --assembly ${hic}_10kb.assemblies.txt  --no-clustering  --prob 0.95
```
---

Then we implement ICE correction following [Imakaev, Maxim et al.](https://www.nature.com/articles/nmeth.2148) and this step has beed packaged in HiC-Pro with one-line command as `HiC-Pro -i raw-matrix -o ice-matrix -c config-hicpro.txt -s ice-norm`.

After the corrections by NeoLoopFinder and ICE, we then condense the chromatin information in Hi-C data.
The defualt way for Hi-C dimention reduction is Singular Value Decomposition (SVD).

```
# used to reduce Hi-C dimention
hic_mat = get_hic_mat(
					data_dir=data_dir,
					drop_rate=hic_drop_rate, #default=0.0
					reduce=hic_reduce, #default='svd'
                            	reduce_dim=hic_reduce_dim, #default=5
                            	resolution=resolution, #default='10Kb', depends on your own Hi-C data
                            	type=hic_type, #default='ice'
                            	norm=hic_norm #default='log')
                            		
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

Now we get the reduced Hi-C data as below:

| gene_name       | HiC-1       | HiC-2       | HiC-3       | HiC-4       | HiC-5       |
|-----------------|-------------|-------------|-------------|-------------|-------------|
| OR4F5           | 1           | 0.000160124 | 0.224521596 | 0.240359624 | 0.48014354  |
| SAMD11          | 0.983033738 | 0.022808685 | 0.215826884 | 0.232252798 | 0.483052779 |
| NOC2L           | 0.982147132 | 0.037410285 | 0.187267998 | 0.224593303 | 0.464582689 |
| KLHL17          | 1           | 0.000160124 | 0.224521596 | 0.240359624 | 0.48014354  |
| PLEKHN1         | 1           | 0.000160124 | 0.224521596 | 0.240359624 | 0.48014354  |
| ISG15           | 0.990748546 | 0.013588959 | 0.217608283 | 0.228825281 | 0.481511407 |
| AGRN            | 0.974279543 | 0.055531025 | 0.190370936 | 0.218615449 | 0.485405181 |
| C1orf159        | 0.96152902  | 0.062132207 | 0.185860829 | 0.212024252 | 0.454002735 |
| TTLL10          | 0.991588209 | 0.010710517 | 0.215210022 | 0.23584837  | 0.479124017 |
| SDF4            | 0.974902568 | 0.039697561 | 0.176009104 | 0.226073947 | 0.49244549  |

#### 1.2 Other omics data

 ```
 code for other omics-data processing 
 ```

#### 1.3 PPI graph construction

Then,we read the PPI data and transform it into a graph through the following commands.

 ```
 # Load PPI data
 ppi_mat = get_ppi_mat(ppi, drop_rate=ppi_drop_rate, from_list=False, random_seed=random_seed, pan=pan) if ppi else None
 edge_index, edge_dim, edge_feat = construct_edge(ppi_mat)
 ```
 
Above functions are shown as below:

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
The basic properties of this graph (on MCF7 ) will be:

- number of edges: 273,765
- number of nodes: 16,165
- feature of any node (e.g., BRCA1): 

| gene_name       | ATAC       | CTCF-1      | CTCF-2      | CTCF-3      | H3K4me3-1   | H3K4me3-2    | H3K27ac-1    | H3K27ac-2    | Means-SNV    | Means-CNV     | Hi-C-1     | Hi-C-2     | Hi-C-3    | Hi-C-4     | Hi-C-5  |      
|------------------|------------|-------------|-------------|-------------|--------------- |----------------|---------------|----------------|----------------|------------------|-----------|------------|-----------|-----------|----------|
| BRCA1              | 0.8721      | 0.4091       | 0.2511        | 0.3964       | 0.8850           |0.9591             |0.9387           |0.8595             |0.6185            | 0.2460               |0.8972     |0.1935       |0.0946      |0.1747     |0.5598     |


Now the input data is prepared. Let'a go for the model training!

 
## 2. Model training

This section demonstrates the GAT-based model architecture of CGMega and how to train CGMega.

---

### CGMega framework

<div align="center"><img width="50%" src="https://github.com/NBStarry/CGMega/blob/main/img/model_architecture.png"/></div>

---

### Hyperparameters

- To reduce the number of parameters and make training feasible within time and resource constraints, the input graphs were sampled using neighbor sampler. The subgraphs included all first and second order neighbors for each node and training was performed on these subgraphs.
- The learning rate is increased linearly from 0 to 0.005 for the first 20% of the total iterations.
- warm-up strategy for learning rate is employed during the initial training phase.
- To prevent overfitting and over smoothing, an early stop strategy is adopted. If the model's performance on the validation set dose not improve for a consecutive 100 epochs, the training process stops.
- Dropout is used and the dropout rate is set to 0.1 for the attention mechanism and 0.4 for the other modules.
- Max pooling step size is 2. After pooling, the representation had 32 dimensions.

---

### System and Computing resources 

| Item     |   Details   |    
|------------------|-----------------------|
| System             | Linux or Windows    |
| RAM  | 256G                       |
| GPU   | NVIDIA GeForce RTX 3090 (24G) |
| Time  | 0.5~1h                    |

 ```note
The above table reports our computing details during CGMega development and IS NOT our computing requirement.

If your computer does not satisfy the above, you may try to lower down the memory used during model training by reduce the sampling parameters, the batch size or so on. 

We are going to test CGMega under more scenarios like with different models of GPU or memory limits to update this table.

 ```

## 3. Interpretation

After prediction, you can do analyses as following to interpret your results.

#### 1. identify the gene module

For each gene, GNNExplainer identified a subgraph G that is most influential to the preiction of its identity from both topological integration and multi-omic information.
This subgraph consists of two parts: 
- i) the core subgraph that consists of the most influential pairwise relationships for cancer gene prediction, and
- ii) the 15-dimension importance scores that quantify the contributions of each gene feature to cancer gene prediction.

```note
The above module identification is calculated at the level of individual genes.
High-order gene modules reported in our work are constructed with the pairwise-connected individual gene modules.
```

#### 2. calculate the Representative Feature

According to the feature importance scores calculated by GNNExplainer, we defined the representative features (RFs) for each gene as its features with relatively prominent importance scores.
In detail, for a given gene, among its features from ATAC, CTCF, H3K4me3 and H3K27ac, SNVs, CNVs and Hi-C, if a feature is assigned with importance score as 10 times higher than the minimum score,it will be referred to as the RF of this gene.
A graphic illustration is shown as below:

<div align="center"><img width="50%" src="https://github.com/NBStarry/CGMega/blob/main/img/RF_calculation.png"/></div>

#### 3. explore the relationships between different gene modules

CGMega serves as a tool to help researchers explore the relationships between individual modules of genes. Such kind of view of high-order gene modules may also helps to find something new.
This is also useful especially when we do not know how to integrate some isolated knowledges into a whole even they are already well-studied under different cases. 

For example, BRCA1 and BRCA2 these two genes act as different roles in common pathway of genome protection and this also exhibited on the topology of their gene modules.
In brief, BRCA1, as a pleiotropic DDR protein working in broad stages of DNA damage response (DDR), was also widely connected with another 20 genes. 
In contrast, BRCA2, as the mediator of the core mechanism of homologous recombination (HR), was connected with other genes via ROCK2, an important gene that directly mediates HR repair.
Moreover, SNV was the RF for both BRCA1 and BRCA2. We also observed a high-order gene module combined from BRCA1 gene module and BRCA2 gene module through three shared genes including TP53, SMAD3 and XPO1.

<div align="center"><img width="50%" src="https://github.com/NBStarry/CGMega/blob/main/img/example.png"/></div>

source: `{{ page.path }}`
