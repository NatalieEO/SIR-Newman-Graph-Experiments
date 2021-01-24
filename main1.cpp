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
	g.ReadGraph("wow.gph");
	Graph g2;
	g2.formGiantComponent(g, 5000);
	g2.recordGraph("copy.gph");
}