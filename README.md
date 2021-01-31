# Summary
This project explores the effect of triangle forming relationships vs clustering in graphs. 

## Sgraph.hpp
This file defines the Graph class and contains a variety of functions to create different types of graphs, save them to a text file, and run SIR/SIR-like experiments on said graph.

### Create a graph 

- void Graph::createGraph(Vsize numberNodes, int gSize, int avgNodeDegree)
    - Not a good function name or function...
    - If the group size is 0, then a random graph is generated
    - Otherwise, this creates a triangle-free graph subdivided into groups.
        - Each group is divided into two halves and edges are created in bipartite fashion within the group.
        - To fulfill the average node degree, edges are then randomly added to nodes outside the group
        - Third argument should be > half the group size
        
- void Graph::TriangleFreeGraph(Vsize numberNodes, int gSize, int inDeg, int avgNodeDegree) {
    - Similar to createGraph but only creates triangle-free graphs and the number of edges within the group is controlled by inDeg.
    - Because a triangle-free graph is achieved in a bipartite manner, the inDeg can't be greater than half the group size.
    
- void TriangleFreeGraph(Vsize numberNodes, int gSize, int inDeg, int avgNodeDegree, int lastEdgeOut)
    - Similiar to TriangleFreeGraph, but only for when one edge is going out of the group per node. ie AvgNodeDegree = inDeg + 1 
    - Last input 'lastEdgeOut' is the probability, out of 10000, that the node's 1 edge going out the group will form
    
- void Graph::Newman(Vertex numNodes, int avgDeg, double cluster, int groupSize)
    - Creates Newman style graph
    
- void Graph::BipNewman(Vertex numNodes, int avgDeg, double cluster, int groupSize)
    - Creates a triangle-free Newman graph by splitting groups in half and adding edges in bipartite manner.
    
### Run experiment on graph
- Giant component approach
    - void formGiantComponent(Graph &g, int p)
        - With input graph g, create another graph by going through each edge and with probability p (out of 10,000), keep an edge.
    - long int Graph::giantComponent()
        - Once giant component has been formed, call giantComponent() to return the size of the largest component found within the graph.
        
- SIR approach
    - long int Graph::SIR(Vertex Seed[], int size, int p)
        - First argument is set of initially infected nodes
        - Second argument is size of set
        - Third argument is probability p out of 10,000 of infecting another node
    - long int Graph::SIRv2(Vertex Seed[], int size, int chance)
        - Same as SIR, but nodes will stay infectious for multiple days. Number of days infected varies for each node.
        
### How to run trials
Note that graph files have .gph extension and trial results are csv files
1. Create and record a graph using one of the following files:
    - genNewmanGraph.cpp
        - Command line arguments: [file path to store graph file to, size of graph, average degree of node, clustering coefficient (double), group size]
    - genBipNewmanGraph.cpp
        - Command line arguments: [file path to store graph file to, size of graph, average degree of node, clustering coefficient (double), group size]
    - genTriangleFreeGraph.cpp
        - Command line arguments: [file path to store graph file to, size of graph, group size, average degree of node within its group, overall average degree of node]
    - genTriangleGraph.cpp
        - Command line arguments: [file path to store graph file to, size of graph, group size, average degree of node within its group, overall average degree of node]
2. Using recorded graph created from step one, use one of the following files to run several trials and save results.
    - runSIRexp.cpp
        - Command line arguments: [file path to graph, file path to save trial results to, group size, average degree of nodes in graph]
        - This was typically ran on a graph with 1 million nodes
        - Runs SIR experiment with infection rates between 100/10000 and 10000/10000, incrementing by 100
        - Each infection rate will be attempted until there are 30 successful runs or until 1000 attempts have been made. Which ever comes first
            - A successful run is when the number of infected nodes is greater than 1000 at the end of the experiment
    - runGiantComponentExp.cpp
        - Command line arguments: [file path to graph, file path to save trial results to, group size, average degree of nodes in graph]
        - This was typically ran on a graph with 1 million nodes
        - Runs giant component experiments with the probability of keeping an edge between 100/10000 and 10000/10000, incrementing by 100
        - Each probability is done 10 times
