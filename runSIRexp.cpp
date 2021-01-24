#include "SGraph.hpp"
#include "Set.hpp"
#include <string>
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
	file << "NumNodes,GroupSize,AvgDeg,NumSeeds,InfRate,NumInfected";
	int gs = atoi(argv[3]);
	int d = atoi(argv[4]);

	Vertex Seed[1] = {0};
	long int numInfected;
	int successRuns = 0;
	int totalAttempts = 0;
	for(int i = 0; i <= 10000; i = i + 100){
		while(successRuns < 30 && totalAttempts < 1000){
			Seed[0] = rand() % g.Size();
			numInfected = g.SIR(Seed, 1, i);
			if(numInfected >= 1000){
				file << '\n' << g.Size() << ", " << gs  << ", " << d << ", " << 1 <<", " << i << ", " << numInfected;
				successRuns++;
			}
			totalAttempts++;
		}
		totalAttempts = 0;
		successRuns = 0;
	}
}