#include "SGraph.hpp"
#include "Set.hpp"
#include <string>
#include <limits>
#include <iostream>
#include <fstream>

int main(int argc, char** argv){
    using namespace std;
    srand(time(NULL));
    
    Vsize size = atoi(argv[2]);
    int groupSize = atoi(argv[3]);
    int inDegree = atoi(argv[4]);
    int avgDegree = atoi(argv[5]);
    
	Graph g;
	g.TriangleFreeGraph(size, groupSize, inDegree, avgDegree);
	g.recordGraph(argv[1]);
}