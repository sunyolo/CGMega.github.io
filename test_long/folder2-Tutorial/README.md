# Tutorial

This tutorial demonstrates how to use CGMega functions with a demo dataset. Once you are familiar with CGMega’s workflow, please replace the demo data with your own data to begin your analysis.

#What the tutorial covers

## 1. How to orepare input data

We recommend getting started with CGMega using the provided demo dataset. When you want to apply CGMega to your own multi-omics dataset, please refer to the following tutorials to learn how to prepare input data. Overall, the input data consists of two parts: the graph, constructed from PPI and the node feature including compressed Hi-C features, SNVs, CNVs frequencies and epigenetic densities.

### PPI graph construction

### Hi-C data embedding

### Other omics data

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
