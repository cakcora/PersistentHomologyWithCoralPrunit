# Persistent Homology with CoralTDA and PrunIt

We develop two new, remarkably simple but effective algorithms to compute the exact persistence diagrams of large graphs to address the computational complexity constraints in TDA. First, we prove that $(k+1)$-core of a graph $G$ is sufficient to compute its $k$ th persistence diagram, $PD_k(G)$. Second, we introduce a pruning algorithm for graphs to compute their persistence diagrams by removing the dominated vertices. Our experiments on large networks indicate that the new approach can achieve computational gains up to 95%.

This notebook shows the PrunIt in action. We find the dominated vertices and delete them from the graph. Our activation function is node degree, and we are using superlevel filtration. Prunit assumes that the dominating vertices are activated before the dominated vertices. Degree in superlevel filtration satisfies this condition. 

Our experiments were coded in R, which is not optimized for loops. We are giving the Python code below, which is easier to read and follow. Sub/super level filtration is not implemented in python yet, but you can call this linked code snippet (to be linked) to send a graph to R and use the TDA library to get its persistence diagram.

```python
import networkx as nx
import time
import random

```


```python
# finds the dominated vertices in a graph. We use degree as the default activation function
# and use sub level filtration. The function does not alter the graph - merely returns a list of dominated
# vertices
def find_dominated_vertices(graph):
    dominated = set()
    for k in graph:
        vec_1_neighbors = set(graph.neighbors(k))
        vec_1_neighbors.add(k)
        len1 = len(vec_1_neighbors)
        for j in graph.neighbors(k):
            vec_2_neighbors = set(graph.neighbors(j))
            vec_2_neighbors.add(j)
            len2 = len(vec_2_neighbors)
            intersection = vec_1_neighbors.intersection(vec_2_neighbors)
            len3 = len(intersection)
            if len3 == len1 and len3 == len2:
                if k not in dominated:
                    dominated.add(j)
            else:
                set_diff1 = vec_1_neighbors.difference(vec_2_neighbors)
                len4 = len(set_diff1)
                if len4 == 0:
                    dominated.add(k)
                else:
                    set_diff2 = vec_2_neighbors.difference(vec_1_neighbors)
                    len5 = len(set_diff2)
                    if len5 == 0:
                        dominated.add(j)
                        
    return dominated
```


```python
# mimics the filtration in Persistent Homology and returns a new graph that
# is computed as the vertices whose activation values (i.e., degrees) are
# higher than a given degree threshold parameter.
def createFiltrationGraph(graph_reduced,graph_full,threshold):
    remove = [node for node,degree in dict(graph_full.degree()).items() if degree <threshold]
    g2 = graph_reduced.copy()
    g2.remove_nodes_from(remove)
    return g2
```


```python
random.seed(10)
# create a graph in three ways
graph = nx.Graph()
# manually input edges
edge_list = [(1,2), (1,3), (1,4), (2,3), (2,5), (3,4), (3,5), (4,6), (5,6)]
graph.add_edges_from(edge_list)
# use an erdos-renyi graph
graph= nx.erdos_renyi_graph(1000,0.001)
# or read a graph from the local disk.
#graph=nx.read_edgelist("Email-Enron.txt")
```


```python
# we will use an interval for filtration. As we are using degree,
# we will use min degree and max degree of the graph as two end points
d = [val for (node, val) in graph.degree()]
maxD = max(d)
minD = min(d)

# degrees in the network may run up to 1000. We use a step size to reduce the number of filtrations
# note that increasing step size has no impact other than extending the run time.
if (maxD-minD)>100:
    step = int((maxD-minD)/100)
else:
    step =1

print('Step size:', step)    
print("Original graph:",graph.number_of_nodes(), graph.number_of_edges())

```

    Step size: 1
    Original graph: 1000 502
    


```python
# our method, PrunIt, eliminated dominated vertices in the very beginning
# unlike strong collapse, PrunIt does not need to check domination at every filtration step.
# as a result, PrunIt is faster than Strong Collapse even when there is no domination.

# note that eliminating vertices at the beginning may shorten the filtration interval as well,
# because some vertices are eliminated and the degrees of the remaining vertices decrease. As a result
# the max degree changes.

######PrunIT - first step: finding and eliminating dominated vertices##########
t3 = time.process_time()
dominated_vertices = find_dominated_vertices(graph)
# we must retain degree information in the original graph
graph_full=graph.copy()
graph.remove_nodes_from((dominated_vertices))
t4 = time.process_time()
print(len(dominated_vertices)," dominated vertices")
```

    305  dominated vertices
    


```python
######PrunIT - second step: using the updated graph in filtration##########
print("Graph nodes and edges:",graph.number_of_nodes(), graph.number_of_edges())
d = [val for (node, val) in graph.degree()]
maxD = max(d)
minD = min(d)
print("Prunit filtration range: ",maxD,minD)
# filtration starts
t5= time.process_time()
totalSimplexCountPrunIt=0
for i in range(maxD, minD-1, -step):
    g2 = createFiltrationGraph(graph,graph_full,i)
    print(i, "th filtration graph nodes and edges:",g2.number_of_nodes(), g2.number_of_edges())
    simplex0 = g2.number_of_nodes()
    simplex1 = g2.number_of_edges()
    simplex2 = sum(nx.triangles(g2).values())
    totalSimplexCountPrunIt+=(simplex0+simplex1+simplex2)^3
t6 = time.process_time()
```

    Graph nodes and edges: 695 197
    Prunit filtration range:  4 0
    4 th filtration graph nodes and edges: 20 8
    3 th filtration graph nodes and edges: 88 46
    2 th filtration graph nodes and edges: 257 197
    1 th filtration graph nodes and edges: 331 197
    0 th filtration graph nodes and edges: 695 197
    


```python
print("Time costs excluding Persistent Homology computations on simplicial complexes:")
print("\t Prunit takes:", (t4-t3)+(t6-t5),"secs")

print("Persistent Homology cost estimates (number of simplices)^3")
print("\t Prunit will create simplicial complexes with :", totalSimplexCountPrunIt," simplices (<=3)")
```

    Time costs excluding Persistent Homology computations on simplicial complexes:
    	 Prunit takes: 0.03125 secs
    Persistent Homology cost estimates (number of simplices)^3
    	 Prunit will create simplicial complexes with : 2043  simplices (<=3)
    

