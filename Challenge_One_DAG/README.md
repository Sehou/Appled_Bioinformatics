# Node_depth_DAG

Depth First Search (DFS) and Breadth First Search (BFS) are two standard algorithms that can be used to find the depth of a node i.e. the  minimum number of edges connecting the given node with the root node in a directed acyclic graph (DAG). 

Node_depth_DAG  is a python script that finds the depth of a specific node n given the adjacency matrix of a DAG in CSV/txt file and the root node. It verifies if there is only one root in the DAG, if the provided graph is acyclic, and if the provided root node provided matched the one it finds in the graph. 

**WARNING:** The adjacency matrix should strictly  be a (0,1) matrix and each row in the file should be separated by a comma. Otherwise the delimiter should explicitly be specified. The adjacency matrix is defined such that a non-zero element *Aij* indicates an edge from *i* to *j*. The node should be an integer and based on 0 index (as indexing in python start from 0) i.e node 1 will literally be provide as 0.

Some test cases are provided along where the right format can be found. 



## Installation

This script does not require any installation. However, its repository must be cloned before usage. This can be done as follow:

```bash
git clone https://github.com/Sehou/Applied_Bioinformatics.git
```



## Usage

The following command can be run in the command line.

```python
python Node_depth_DAG.py <options>
```

The options -i for the CSV/txt file containing the adjacency matrix, -r for root node and -n for the node of interest are mandatory. All the available options can be see by running: 

```bash
python Node_depth_DAG .py --help
```

## Output

The program print out the output to the command line.