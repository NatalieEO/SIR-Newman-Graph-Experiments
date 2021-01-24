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

// 	Graph g;
// 	g.createTriGraph(10000, 10, 7, 10);
// 	g.recordGraph("halp.gph");

//     Graph g2;
// 	g2.ReadGraph("halp.gph");


//     //Vertex Seed[] = {500000, 600000, 700000};
// 	Vertex Seed[3] = {0, 0, 0};
//     ofstream file;
//     file.open("trying.csv");
// 	file << "NumNodes,GroupSize,AvgDeg,NumSeed,InfRate,NumInfected";
//     int numSeeds = 3;
//     for(int i = 0; i <= 7000; i = i + 100){
//         for(int j = 0; j < 1; j++){
//     		for(int k = 0; k < 3; k++){
//     			do{
//     				Seed[k] = rand() % 10000;
//     			}while(Seed[k] == Seed[(k+1)%3] || Seed[k] == Seed[(k+2)%3]);	
//     		}
// 		    //cout << Seed[0] << " " << Seed[1] << " " << Seed[2];
// 		    Graph g3;
//             g3.formGiantComponent(g2, i);
//             long int numInfected = g3.giantComponent(0);
//             file << '\n' << 1000000 << ", " << 10 << ", " << 6 << ", " << 3 << ", " << i << ", " << numInfected;
//         } 
//     }
    Graph g;
    g.BipNewman(10000, 5, 0.5, 10);
}