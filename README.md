# Summary
This project explores the effect of triangle forming relationships vs clustering in graphs. 

## Sgraph.hpp
This file defines the Graph class and contains a variety of functions to create different types of graphs, save them to a text file, and run SIR/SIR-like experiments on said graph.

### Create a graph functions

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
    
### Run experiment on graph functions
- Giant component approach
    - void formGiantComponent(Graph &g, int p)
        - With input graph g, create another graph by going through each edge and with probability p (out of 10,000), keep an edge.
    - long int Graph::giantComponent()
        - Once giant component has been formed, call giantComponent() to return the size of the largest component found within the graph.
        
- SIR approach
