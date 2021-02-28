#include "SGraph.hpp"
#include "Set.hpp"
#include <string>
#include <stdlib.h>
#include <limits>
#include <iostream>
#include <fstream>

//g++-5 -std=c++14 main.cpp
int main(int argc, char** argv){
    using namespace std;
    srand(time(NULL));


	Graph g;
	g.ReadGraph(argv[1]);
    ofstream file;
    file.open(argv[2]);
    file << "NumNodes,GroupSize,AvgDeg,InfRate,NumInfected";

    int gs = atoi(argv[3]);
    int d = atoi(argv[4]);
    long int numInfected;
    for(int i = 1000; i <= 2000; i = i + 100){
        for(int j = 0; j < 10; j++){
 	        Graph g2;
            g2.formGiantComponent(g, i);
            numInfected = g2.giantComponent();
            file << '\n' << g.Size() << ", " << gs  << ", " << d << ", " << i << ", " << numInfected;
        } 
    }
}