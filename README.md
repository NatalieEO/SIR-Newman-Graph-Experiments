# Summary
This project explores the effect of triangle forming relationships vs clustering in graphs. 

# Setup
To run the Jupyter Notebooks, this is confirmed working with Python 3.7.0 and pip 10.0.1. You will need to install:
- pip install pandas
- pip install matplotlib
- pip install plotly==3.10.0

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
    - Graph is formed by:
            - Creates a graph with number of nodes n.
            - Divides nodes into groups of size gSize.
            - Each group is divided in half. A node in one half will attempt to form an edge to each node in the second half to reach desired inDeg
            - Once all nodes reach desired inDeg, form edges randomly to nodes outside of group to get desired avgNodeDegree. 
    
- void TriangleFreeGraph(Vsize numberNodes, int gSize, int inDeg, int avgNodeDegree, int lastEdgeOut)
    - Similiar to TriangleFreeGraph, but only for when one edge is going out of the group per node. ie AvgNodeDegree = inDeg + 1 
    - Last input 'lastEdgeOut' is the probability, out of 10000, that the node's 1 edge going out the group will form
     - Graph is formed by in same manner as above, but a node will only attempt to form an edge outside of the group once with probability p.
    
- void Graph::createTriGraph(Vsize numberNodes, int gSize, int avgInDegree, int avgNodeDegree)
    - Creates a graph with triangles subdivided into groups.
    - avgInDegree should be lower than the average node degree
    - Graph is formed by:
        - Creates a graph with number of nodes n.
        - Divides nodes into groups of size gSize.
        - Within a group, each node will attempt to form an edge, with probability p (calculated), with another node in the same group to reach avgInDegree.
        - Once all nodes reach desired avgInDegree, form edges randomly to nodes outside of group to get desired avgNodeDegree. 

- void Graph::createTriGraph(Vsize numberNodes, int gSize, int avgInDegree, int avgNodeDegree, int lastEdgeOut)
    - Similar to createTriGraph, but but only for when one edge is going out of the group per node. ie AvgNodeDegree = inDeg + 1
    - Last input 'lastEdgeOut' is the probability, out of 10000, that the node's 1 edge going out the group will form
    - Graph is formed in same manner as above, but a node will only attempt to form an edge outside of the group once with probability p.
    
- void Graph::Newman(Vertex numNodes, int avgDeg, double cluster, int groupSize)
    - Creates Newman style graph:
        - Repeatedly create different groups of size groupSize with randomly chosen nodes.
        - Attempt to join attempt to join an edge with probability p between each node within group once
            - Probability p calculated using clustering coefficient cluster, average node degree avgDeg, and group size groupSize.
        - Keep forming groups and adding edges until needed number of edges to have average node degree of avgDeg is reached.
    
- void Graph::BipNewman(Vertex numNodes, int avgDeg, double cluster, int groupSize)
    - Creates a triangle-free Newman graph using method described above. However, now each group is splity in half and edges are added in bipartite manner.
        - The probability p that each node can only connect with nodes in other half of group is almost double than that of p in regular Newman graph.
    
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
        
## How to run trials
Note that graph files have .gph extension and trial results are csv files
1. Create and record a graph using one of the following files:
    - genNewmanGraph.cpp
        - Command line arguments: [file path to store graph file to (.gph), size of graph, average degree of node, clustering coefficient (double), group size]
    - genBipNewmanGraph.cpp
        - Command line arguments: [file path to store graph file to, size of graph, average degree of node, clustering coefficient (double), group size]
    - genTriangleFreeGraph.cpp
        - Command line arguments: [file path to store graph file to, size of graph, group size, average degree of node within its group, overall average degree of node]
    - genTriangleGraph.cpp
        - Command line arguments: [file path to store graph file to, size of graph, group size, average degree of node within its group, overall average degree of node]
2. Using recorded graph created from step one, use one of the following files to run several trials and save results. Check 1milData description below to see naming format of .csv files
    - runSIRexp.cpp
        - Command line arguments: [file path to graph, file path to save trial results to (.csv), group size, average degree of nodes in graph]
        - This was typically ran on a graph with 1 million nodes
        - Runs SIR experiment with infection rates between 100/10000 and 10000/10000, incrementing by 100
        - Each infection rate will be attempted until there are 30 successful runs or until 1000 attempts have been made. Which ever comes first
            - A successful run is when the number of infected nodes is greater than 1000 at the end of the experiment
    - runGiantComponentExp.cpp
        - Command line arguments: [file path to graph, file path to save trial results to (.csv), group size, average degree of nodes in graph]
        - This was typically ran on a graph with 1 million nodes
        - Runs giant component experiments with the probability of keeping an edge between 100/10000 and 10000/10000, incrementing by 100
        - Each probability is done 10 times
        
    - The .csv files created from the 2 files above follow
        
### Data
The data generated is kept in the "Trial Data" folder. 
- 1milData 
    - This folder contains all the .csv files generated from the SIR or giant component experiments.
    - Sub folders say size of the graph, number of seeds, and if seeds stay in infectious for multiple days
    - The csv files follow a similar naming format. 
        - Type of graph
        - 1SIR30 -> 1(seed), 30(trials) per infection rate
        - Size of graph
        - GS10 -> Group size = 10
        - AD5 -> Average Degree of a node = 5 
        - OR AD26-15 -> Average degree of a node = 26 but average degree of a node to other nodes within same group = 15
        - C02 -> For Newman style graphs only, Cluster coefficient = 0.2
        
- Big Graph Experiments
    - Contains all the jupyter notebooks created so far.

