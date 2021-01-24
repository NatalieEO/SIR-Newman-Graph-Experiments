/* 
 * File:   SGraph.hpp
 * Author: maleq
 *
 * Created on November 13, 2008, 4:35 PM
 */

#ifndef _GRAPH_HPP
#define	_GRAPH_HPP

#include <map>
#include <iostream>  
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <cmath>
#include <list>
#include <queue>
#include <time.h>

/* 
 * File:   SGraph.hpp
 * Author: maleq
 *
 */
 
 
#include "utility.hpp"
#include "Set.hpp"

#ifndef NULL
	#define  NULL 0
#endif

#define  FNAME_LEN 512	//maximum length of graph file name
#define  MAX_LABEL 6

using namespace std;

typedef  unsigned long long int  ULLI;
typedef  unsigned long int  ULI;
typedef  unsigned long int  Vertex;
typedef  unsigned long int  Vsize;
typedef  unsigned short int Word;
typedef  unsigned char      Byte;
typedef  Vertex*  Vlist;


class Graph {
public:
    Vsize    n;			// # of vertices
	Vsize	 maxd;		// maximum degree
    Set		*nlist;		// neighbor list. nlist[i] is the neighbor list of vertex i
    Byte    *label;		// labels of vertices	
	ULI		*name;		// original ID of the vertices
    char	gfname[FNAME_LEN];	// graph file name: set in ReadGraph
    //int objectCount;

	Graph() {
		nlist = NULL;
		label = NULL;
		name = NULL;
		maxd = 0;
	}
	
	~Graph() {
		FreeGraph();
	}
	
	void BFS(Vertex source, Graph &g);
	int SIR(Vertex Seed[], int size ,Graph &g, int chance);
	void createGraph(Vsize numberNodes, int groupSize, int avgNodeDegree);
	void bipartite(int startNode, int groupSize);
	void randomEdges(int i, int groupSize, float half);
	void createTriGraph(Vsize numberNodes, int gSize, int avgInDegree, int avgNodeDegree);
	
	Vsize ReadGraph(const char *fname, const char *demfile);
	Vsize Deg(Vertex v) {return nlist[v].Size();}
	ULI Name(Vertex v) {return name[v];}
	Byte &Label(Vertex v) {return label[v];}
	Vsize Size() {return n;}
	Vertex &operator()(Vertex v, Vsize i) {return nlist[v][i];}
	Set &operator[](Vertex v) {return nlist[v];}
	void FreeGraph();
	Vsize MaxDeg();
	void Print();
	void Sort();
	void CheckGraph();
	void rNeighborhood(Vertex v, int r, Set &Vr);	// computes the set of nodes in r-neighborhood
	void rNeighborhoodRec(Vertex v, int r, Set &Vr); // recursive func - used by rNeighborhood
	void rNeighborhoodGraph(Vertex v, int r, Graph &g); // Computes r-neighborhood graph
};



Vsize Graph::ReadGraph(const char *fname, const char *demfile=NULL) 
{
	ULI     u, v, temp;					// vertices
	Vsize   i, k, tmpdeg, dummy;
	Vertex  vidx;					// vertex index
	double  wt;
	map<ULI, Vertex> vindex;
	//objectCount = 0;

	strcpy(gfname, fname);
	ifstream ifp(fname);
	if(!ifp.is_open()) {
		cout << "\% Cannot open " << fname << endl;
		exit(1);	
	}

	//cout << "Reading graph file: " << fname << " ... "; cout.flush();
	ifp >> n;									// read first line -- the number of vertex
    nlist = new Set[n];
	name = new ULI[n];
	vidx = 0;
	
 	for (i=0; i<n; i++) {
    	ifp >> u >> tmpdeg;						// read nodes and its degree

		if (!vindex.count(u)) {					// if new vertex, map it 
			name[vidx] = u;
			vindex[u] = vidx++;
		}
		
		nlist[vindex[u]].init(tmpdeg);
	
		for(k=0; k<tmpdeg; k++) {
      		ifp >> v >> wt >> temp;					// read the adjacent nodes

			if (!vindex.count(v)) {
				name[vidx] = v;
				vindex[v] = vidx++;
			}
			
			nlist[vindex[u]].insert(vindex[v]);		// add the edge (u, v) to the graph
   		}
	}
	ifp.close();
	if (demfile==NULL) return n;
   
		// read dem file and create labels for the virtices

	cout << "done. \n\% Reading the demographic file: " << demfile << " ... "; cout.flush();
	
	label = new Byte[n];
	for (i=0; i<n; i++)
    	label[i] = 1;
	
	ifp.open(demfile);
	if(!ifp.is_open()) {
        cout << "\% Cannot open demographic file: " << fname << endl;
        exit(1);	
	}

	double vul;
	//ifp.ignore(1000,'\n');				// ignore first line: instructions
	while(!ifp.eof()) {
		//ifp.ignore(1000,'\t');			// ignore first column
		ifp >> u >> vul;			// read id and age
		ifp.ignore(1000,'\n');			// ignore the rest of the line
		if (vindex.count(u)) {
			label[vindex[u]] = (vul < 0.2) ? 0 : ((vul > 0.8) ? 2 : 1);  // label using ages only
		}
	}
	ifp.close();
	
	//for (i=0; i<n; i++)
    //	if (label[i] == 126) cout << "missing age of person " << name[i] << endl;
	
	cout << "done" << endl;
 	return n; 
}



void Graph::FreeGraph()
{
	if (nlist)
		for (Vsize i=0; i<n; i++) 
			nlist[i].destroy();
	FreeMem(nlist);
	FreeMem(label);
	FreeMem(name);
}



Vsize Graph::MaxDeg() 
{
	if (maxd==0) {
		for (Vertex v=0; v<n; v++)
			if (maxd < nlist[v].Size())
				maxd = nlist[v].Size();
	}
	return maxd;
}



void Graph::Print()
{
	Vertex i;
	Vsize k;
	
	if (nlist==NULL) return;
	cout << "\% Graph: \n" << n << endl;
	for (i=0; i<n; i++) {
		cout << "\% " << name[i] << "  " << nlist[i].Size() << endl; 
		for (k=0; k<nlist[i].Size(); k++)
			cout << "\% \t" << name[nlist[i][k]] << endl;
	}
}



// check if the graph is a consistent simple undirected graph
void Graph::Sort()
{
	Vertex v;

	for (v=0; v<n; v++) 
		nlist[v].sort();
}



void Graph::CheckGraph()
{
	Vertex i, v;
	Vsize k;

	cout << "Begin graph checking ..." << endl;

	for (i=0; i<n; i++) {
		if (!Deg(i))
			cout << "Degree-zero node: " << name[i] << endl;
		else	
			nlist[i].sort();
	}
	for (i=0; i<n; i++) {
		for (k=1; k<Deg(i); k++) {
			if (nlist[i][k-1]==nlist[i][k])
				cout << "Duplicate neighbor: " << name[i] << " " << name[nlist[i][k]] << endl;
		}
		for (k=0; k<Deg(i); k++) {
			v = nlist[i][k];
			if (v > n-1)
				cout << "Unknown neighbor: " << Name(i) << " " << Name(v) << endl;
			else if (v == i)
				cout << "Self-loop: " << Name(i) << " " << Name(v) << endl;
			else if (!nlist[v](i))
				cout << "Oneway neighbor: " << name[i] << " " << name[v] << endl;
		}
		if (label && label[i] > MAX_LABEL - 1) 
			cout << "Large label: " << name[i] << " " << (unsigned int) label[i] << endl;
	}
	cout << "Graph checking done." << endl;
}



// called by rNeighborhood
void Graph::rNeighborhoodRec(Vertex v, int r, Set &Vr)
{
	Vr.set_union(nlist[v]);
	if (r==1) return;
	for (Vsize i=0; i<Deg(v); i++)
		rNeighborhoodRec(nlist[v][i], r-1, Vr);
}


//Find Vr, the set of all nodes, including v, within distance r from v
//g.Sort() must be called before calling this function
void Graph::rNeighborhood(Vertex v, int r, Set &Vr)
{
	if (r<0) return;
	Vr.init(1);
	Vr.insert(v);
	if (r==0) return;
	rNeighborhoodRec(v, r, Vr);
}


//Generates r-neighborhood graph, induced graphs by the nodes in r-neighborhood 
// of vertex v (including node v).
//g.Sort() must be called before calling this function
void Graph::rNeighborhoodGraph(Vertex v, int r, Graph &g) 
{
	Set V;
	
	rNeighborhood(v, r, V);
	g.n = V.Size();		
    g.nlist = new Set[g.n];
	g.name = new ULI[g.n];

	Vsize i;
	for (i=0; i<g.n; i++) {
		g[i].init(Deg(V[i]));
		g[i].intersect_idx(V, nlist[V[i]]);
		g.name[i] = name[V[i]];
	}
}

void Graph::BFS(Vertex source, Graph &g){
	// All vertices are not visited by default
	bool visited[n];
    for(int i = 0; i < n; i++){
        visited[i] = false;
    }
    
    // All vertices have an undetermined distance away from source by default.
    int distanceFromSource[n];
    for(int i = 0; i < n; i++){
        distanceFromSource[i] = 0;
    }
    
    visited[source] = true;
    
    // make another list/queue to keep track of what vertexes to visit from the source
    list<int> queue;
    // the only thing in the queue in the beginning is the source
    queue.push_back(source);
    
    while(!queue.empty()){
    	
    	//Set source equal to what is at front of queue
        source = queue.front();
        int i;
        // Print out the traversal, and distance from source every vertex is.
        cout << "Vertex ID: " << g.Name(source) << " Distance from source: " 
        << distanceFromSource[source] << '\n';
        
        // Delete the first element of the queue
        queue.pop_front();
    	
    	for(i = 0; i < nlist[source].Size(); i++){
    		if(!visited[nlist[source][i]]){
    			distanceFromSource[nlist[source][i]] = distanceFromSource[source] + 1;
    			visited[nlist[source][i]] = true;
    			
    			queue.push_back(nlist[source][i]);
    		}
    	}
	}
}



int Graph::SIR(Vertex Seed[], int size, Graph &g, int chance){
	int* infectionDay = new int[n];
	int numInfected = size;
	int random;
	//Set every node's infection day to -1
	for(int i = 0; i < n; i++){
		infectionDay[i] = -1;
	}
	
	queue<int> Q;
	
	//Put initially infected nodes into the queue
	for(int i = 0; i < size; i++){
		infectionDay[Seed[i]] = 0;
		Q.push(Seed[i]);
	}
	
	Vertex source;
	
	while(!Q.empty()){
		source = Q.front();
		Q.pop();
		nlist[source].print();
		
		for(int i = 0; i < nlist[source].Size(); i++){
			//if node hasn't been infected and the generated chance is less the the probability of getting infected
			//then add it to the queue of infected nodes. Set infection day as 1 plus its parent.
			random = rand() % 10000;
			cout << infectionDay[nlist[source][i]] << endl;
			if(infectionDay[nlist[source][i]] == -1 && (random) < chance){
				infectionDay[nlist[source][i]] = infectionDay[source] + 1;
				Q.push(nlist[source][i]);
				numInfected++;
			}
		}
	}
	return numInfected;
}


void Graph::createTriGraph(Vsize numberNodes, int gSize, int avgInDegree, int avgNodeDegree){
	n = numberNodes;
	nlist = new Set[n];
	name = new ULI[n];
	float probability = (float(avgInDegree)) / (gSize - 1);
	int p = probability * 10000;
	int numberOfGroups = n / gSize;
	
	// initialize size for each node in the graph
	for(Vertex i = 0; i < n; i++){
	    name[i] = i;
	    nlist[i].init(avgNodeDegree * 2);
	}
	// For each group, attempt to join an edge between each vertex and every other 
	// vertex in the group.
	for(int i = 0; i < numberOfGroups; i++){
		for(int j = 0; j < gSize-1; j++){
			for(int k = j+1; k < gSize; k++){
				if(rand() % 10000 < p){
					nlist[(i*gSize) + j].Dinsert((i*gSize) + k);
					nlist[(i*gSize) + k].Dinsert((i*gSize) + j);
				}
			}
		}
	}

	// If the last group has fewer nodes than the rest
	if(gSize*numberOfGroups != n){
		Vertex remaining = n - (gSize*numberOfGroups);
		for(int j = 0; j < remaining; j++){
			for(int k = j+1; k < remaining; k++){
				if(rand() % 10000 < p){
					nlist[j+(gSize*numberOfGroups)].Dinsert(k+(gSize*numberOfGroups));
					nlist[k+(gSize*numberOfGroups)].Dinsert(j+(gSize*numberOfGroups));
				}
			}
		}
	}
	
	// To reach the desired overall average node degree, go through each node again
	// and connect edges randomly.
	float half = (avgNodeDegree - avgInDegree) / 2.0;
	for(int i = 0; i < n; i++){
		randomEdges(i, gSize, half);
	}
	
	float sum = 0;
	for (int i = 0; i < n; i++){
		sum = nlist[i].Size() + sum;
	}
	cout << sum / 35.0 << '\n';
}


void Graph::createGraph(Vsize numberNodes, int gSize, int avgNodeDegree){
	
	n = numberNodes;
	nlist = new Set[n];
	name = new ULI[n];
	
	if(gSize == 0){
		//Set up the size of all the neighbor lists
		for(int i = 0; i < n; i++){
	    	name[i] = i;
	    	nlist[i].init(avgNodeDegree * 2);
		}
		
		for(int i = 0; i < n; i++){
			randomEdges(i, gSize, avgNodeDegree/2.0);
		}
	}
	
	else {
		int numberOfGroups = n / gSize;
		int iterations = numberOfGroups * gSize;
		float qi = avgNodeDegree - (gSize / 2.0);
		int setSize = (qi * 3)/2;
		float half = qi / 2;
		//srand(time(NULL));
		
		//Set up the size of all the neighbor lists
		for(int i = 0; i < n; i++){
	    	name[i] = i;
	    	nlist[i].init(setSize);
		}
		
		for(int i = 0; i < numberOfGroups; i++){
			bipartite(i*gSize, gSize);
		}
		if(iterations != n){
			bipartite(iterations, n - iterations);
		}
		
		for(int i = 0; i < n; i++){
			randomEdges(i, gSize, half);
		}		
	}
	
	// float sum = 0;
	// for (int i = 0; i < n; i++){
	// 	sum = nlist[i].Size() + sum;
	// }
	// cout << sum / 100.0 << '\n';
}


void Graph::bipartite(int startNode, int groupSize){
	int midNode = startNode + (groupSize/2);
	int lastNode = startNode + groupSize;
	for(int i = startNode; i < midNode; i++){
		for(int k = midNode; k < lastNode; k++){
			nlist[i].Dinsert(k);
			//Connect node in second half to node from first half.
			nlist[k].Dinsert(i);
		}
	}
}


void Graph::randomEdges(int i, int groupSize, float half){
	int randomNode;     			//The node to add.
	bool prevConnected;
	bool sameGroup;
	bool itself;
	if(groupSize != 0){
		int groupNum = i/groupSize;		//What group number this node belongs to.
		for(int j = 0; j < half; j++){
			randomNode = rand() % n ;
			
			// Keep generating a random node to attach to if...
			// it's been connected to previously
			prevConnected = nlist[i].member(randomNode);
			// Lies in the same group
			sameGroup = (randomNode >= (groupNum * groupSize)) && (randomNode < ((groupNum + 1) * groupSize));
			// Connects to itself
			itself = i == randomNode;
			while(prevConnected || sameGroup || itself){
				randomNode = rand() % n;
				prevConnected = nlist[i].member(randomNode);
				sameGroup = (randomNode >= (groupNum * groupSize)) && (randomNode < ((groupNum + 1) * groupSize));
				itself = i == randomNode;
			}
			
			//Insert edge between the current node and random node selected
			nlist[i].Dinsert(randomNode);
			nlist[randomNode].Dinsert(i);
		}	
	}

	// There are no group subdivisions
	else{
		for(int j = 0; j < half; j++){
			randomNode = rand() % n ;
			prevConnected = nlist[i].member(randomNode);
			itself = i == randomNode;
			
			while(prevConnected || itself){
				randomNode = rand() % n;
				prevConnected = nlist[i].member(randomNode);
				itself = i == randomNode;
			}
			
			nlist[i].Dinsert(randomNode);
			nlist[randomNode].Dinsert(i);
		}
	}
}
#endif	/* _GRAPH_HPP */