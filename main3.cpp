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
	g.TriangleFreeGraph(10000, 20, 9, 13);
	g.recordGraph("TriangleFreeGraph.gph");
}