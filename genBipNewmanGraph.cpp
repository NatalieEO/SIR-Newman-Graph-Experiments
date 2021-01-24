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
    int avgDegree = atoi(argv[3]);
    double clusterC = atof(argv[4]);
    int groupSize = atoi(argv[5]);
    
	Graph g;
	g.BipNewman(size, avgDegree, clusterC, groupSize);
	g.recordGraph(argv[1]);
}