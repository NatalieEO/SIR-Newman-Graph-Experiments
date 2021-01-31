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
	int BFSmod(Vertex source, bool visited[]);
	long int giantComponent();
	long int linearSearch(long int index, bool visited[]);
	long int SIR(Vertex Seed[], int size, int p);
	long int dSIR(Vertex Seed[], int size, double p);
	long int SIRv2(Vertex Seed[], int size, int chance);
	void formGiantComponent(Graph &g, int p);
	void recordGraph(const char *fname);
	void createGraph(Vsize numberNodes, int groupSize, int avgNodeDegree);
	void TriangleFreeGraph(Vsize numberNodes, int gSize, int inDeg ,int avgNodeDegree);
	void TriangleFreeGraph(Vsize numberNodes, int gSize, int inDeg, int avgNodeDegree, int lastEdgeOut);
	void Newman(Vertex numNodes, int avgDeg, double cCoef, int groupSize);
	void BipNewman(Vertex numNodes, int avgDeg, double cCoef, int groupSize);
	bool searchArray(Vertex k, Vertex array[], int size);
	void bipartite(Vertex startNode, int groupSize);
	void bipartiteP(Vertex startNode, int groupSize, int p);
	void randomEdges(Vertex i, int groupSize, int half);
	void createTriGraph(Vsize numberNodes, int gSize, int avgInDegree, int avgNodeDegree);
	void createTriGraph(Vsize numberNodes, int gSize, int avgInDegree, int avgNodeDegree, int lastEdgeOut);
	void gnr(Vsize numberNodes, int avgDegree);
	
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
	long int numEdges = 0;
	for(Vertex i = 0; i < n; i++){
		numEdges = numEdges + nlist[i].Size();
	}
	cout << numEdges << endl;
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


// Write to graph to a separate file
void Graph::recordGraph(const char *fname){
	Vertex i;
	Vsize k;
	
	if(nlist==NULL) return;
	ofstream file;
	file.open(fname);
	file << n << endl;
	for(i = 0; i < n; i++){
		file << name[i] << "  " << nlist[i].Size() << endl;
		for(k = 0; k < nlist[i].Size(); k++){
			file << '\t' << name[nlist[i][k]] << "\t0\t0" << endl;
		}
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
    
    // make another list/queue to keep track of what vertices to visit from the source
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


// Use BFS search to return size of component the source node belongs to and update which nodes have been visited
int Graph::BFSmod(Vertex source, bool visited[]){
    int size = 1;
    visited[source] = true;
    
    // make another list/queue to keep track of what vertices to visit from the source
    list<int> queue;
    // the only thing in the queue in the beginning is the source
    queue.push_back(source);
    int i;
    while(!queue.empty()){
    	
    	//Set source equal to what is at front of queue
        source = queue.front();
        
        // Delete the first element of the queue
        queue.pop_front();
    	
    	for(i = 0; i < nlist[source].Size(); i++){
    		if(!visited[nlist[source][i]]){
    			visited[nlist[source][i]] = true;
    			queue.push_back(nlist[source][i]);
    			size++;
    		}
    	}
	}
	return size;
}


// return node that hasn't been visited yet during BFS
long int Graph::linearSearch(long int index, bool visited[]){
	for(Vertex i = index; i < n; i++){
		if (visited[i] == false){
			return i;
		}
	}
	// If all nodes have been visited already, return -1
	return -1;
}


// Return the size of the giant component
long int Graph::giantComponent(){
    bool visited[n];					// Keeps track whether a node has been visited
    bool searchComplete = false;
    long int index = 0;
    int giantComponentSize = 0;
    int compare;
    
    // All nodes have not been visited in the begining
    for(Vertex i = 0; i < n; i++){
    	visited[i] = false;
    }
    
    visited[0] = true;
    
    while(!searchComplete){
    	// BFSmod returns size of component node is a part of and continues to update the visited array
    	compare = BFSmod(index, visited);
    	if(giantComponentSize < compare){
    		giantComponentSize = compare;
    	}
    	
    	index = index + 1;
    	// Set index to the next unvisited node. 
    	index = linearSearch(index, visited);
    	// if index = -1, then all nodes have been visited already. This will terminate while loop.
    	if(index == -1){
    		searchComplete = true;
    	}
    }
    return giantComponentSize;
}


// SIR simulation. Nodes will stay infectious for only 1 day.
long int Graph::SIR(Vertex Seed[], int size, int p) {
	int infectionDay[n];
	long int numInfected = size;
	int random;
	//Set every node's infection day to -1
	for(Vertex i = 0; i < n; i++){
		infectionDay[i] = -1;
	}
	
	queue<Vertex> Q;
	// Make sure Q is empty before adding any elements
	if(!Q.empty()){
		return -1;
	}
	
	//Put initially infected nodes into the queue
	for(Vertex i = 0; i < size; i++){
		infectionDay[Seed[i]] = 0; 
		Q.push(Seed[i]);
	}
	
	Vertex source;
	int counter;
	while(!Q.empty()){
		source = Q.front();
		Q.pop();
		
		for(int i = 0; i < nlist[source].Size(); i++){
			//if node hasn't been infected and the generated chance is less the the probability of getting infected
			//then add it to the queue of infected nodes. Set infection day as 1 plus its parent.
			random = rand() % 10000;
			if(infectionDay[nlist[source][i]] == -1 && random < p){
				infectionDay[nlist[source][i]] = infectionDay[source] + 1;
				Q.push(nlist[source][i]);
				numInfected++;
			}
		}
	}
	return numInfected;
}

// Same as SIR but p is double type instead of int type.
long int Graph::dSIR(Vertex Seed[], int size, double p){
	int infectionDay[n];
	long int numInfected = size;
	double random;
	//Set every node's infection day to -1
	for(Vertex i = 0; i < n; i++){
		infectionDay[i] = -1;
	}
	
	queue<Vertex> Q;
	// Make sure Q is empty before adding any elements
	if(!Q.empty()){
		return -1;
	}
	
	//Put initially infected nodes into the queue
	for(Vertex i = 0; i < size; i++){
		infectionDay[Seed[i]] = 0;
		Q.push(Seed[i]);
	}
	
	Vertex source;
	int counter;
	while(!Q.empty()){
		source = Q.front();
		Q.pop();
		// nlist[source].print();
		
		for(int i = 0; i < nlist[source].Size(); i++){
			//if node hasn't been infected and the generated chance is less the the probability of getting infected
			//then add it to the queue of infected nodes. Set infection day as 1 plus its parent.
			
			random = Uniform(0.0, 1.0);
			// cout << infectionDay[nlist[source][i]] << endl;
			if(infectionDay[nlist[source][i]] == -1 && random < p){
				infectionDay[nlist[source][i]] = infectionDay[source] + 1;
				Q.push(nlist[source][i]);
				numInfected++;
			}
		}
	}
	return numInfected;
}


// SIR version where nodes will stay infectious for multiple days.
long int Graph::SIRv2(Vertex Seed[], int size, int chance){
	long int infectionDay[n];
	int infDaysLeft[n];					// The number of days of left for a node to stay infectious
	int setInfDays[3] = {3, 4, 5};		// Used to determine the number of days a node will be infectious
	long int numInfected = size;
	int random;							// random = rand() % 10000, if random < chance, node becomes infected.
	int dayP;							// used for determining how many days a node will stay infectious
	//Set every node's infection day to -1
	for(Vertex i = 0; i < n; i++){
		infectionDay[i] = -1;
	}
	
	queue<Vertex> Q;
	// Make sure Q is empty before adding any elements
	if(!Q.empty()){
		return -1;
	}
	
	//Put initially infected nodes into the queue, and set their infection day to zero
	for(Vertex i = 0; i < size; i++){
		infectionDay[Seed[i]] = 0;
		dayP = rand() % 3;
		//Set up the number of days the seed nodes will stay infectious
		infDaysLeft[Seed[i]] = setInfDays[dayP];
		Q.push(Seed[i]);
		//cout << "Vertex " << Seed[i] << " will stay infections for " << infDaysLeft[Seed[i]] << " days" << endl;
	}
	
	Vertex source;
	int counter;
	while(!Q.empty()){
		source = Q.front();
		//cout << "Infected Vertex " <<source << " has " << infDaysLeft[source] << " days left" << endl;
		infDaysLeft[source] = infDaysLeft[source] - 1;
		Q.pop();
	
		for(int i = 0; i < nlist[source].Size(); i++){
			//if node hasn't been infected and the generated chance is less the the probability of getting infected
			//then add it to the queue of infected nodes. Set infection day as 1 plus its parent.
			random = rand() % 10000;
			// cout << infectionDay[nlist[source][i]] << endl;
			if(infectionDay[nlist[source][i]] == -1 && (random) < chance){
				infectionDay[nlist[source][i]] = infectionDay[source] + 1;
				Q.push(nlist[source][i]);
				
				dayP = rand() % 3;
				infDaysLeft[nlist[source][i]] = setInfDays[dayP];
				numInfected++;
			}
		}
		
		if(infDaysLeft[source] != 0){
			Q.push(source);
		}
	}
	return numInfected;
}



// Create another graph using input graph g and deleting edges from it. More flexible to create copy of graph in case you want to
// to do something to the original graph g during the same main call.
void Graph::formGiantComponent(Graph &g, int p){
	n = g.Size();
	nlist = new Set[n];
	name = new ULI[n];
	Set *hasConnection = new Set[n];
	int chance;				// variable to hold rand() % 10000
	
	//copy original graph g
	for(Vertex i = 0; i < n; i++){
		Vsize size = g.nlist[i].Size();
		nlist[i].init(size);
		name[i] = g.name[i];
		for(Vsize j = 0; j < size; j++){
			nlist[i].insert(g.nlist[i][j]);
		}
	}
	
	for(Vertex i = 0; i < n; i++){
		for(Vsize j = 0; j < nlist[i].Size(); j++){
			chance = rand()%10000;
			// Don't attempt to delete the same edge twice. Cannot use name[i] when input graph is created by reading from a file in main
			if(i < nlist[i][j]){
				// Keep an edge with probability p aka Delete edge if rand % 10000 >= p.
				if(chance >= p){
					// Delete edge. Set neighbor that's being deleted to the last neighbor in the list and decrease neighbor list size
					// Decrement j so newly switched neighbor is not skipped in next iteration.
					Vertex temp = nlist[i][j];
					nlist[i][j] = nlist[i][nlist[i].Size() - 1];		
					nlist[i].decSize();
					j--;	
					
					// Find the corresponding neighbor list and switch last neighbor entry with neighbor that's to be deleted.
					// ie nlist[1][2] = 3, find nlist[3][1] and delete.
					Vsize k = 0;
					while(nlist[temp][k] != i){
						k++;
					}
					nlist[temp][k] = nlist[temp][nlist[temp].Size() - 1];
					nlist[temp].decSize();
				}
			}
		}
	}
}


// Create a graph with triangles subdivided into groups. 
void Graph::createTriGraph(Vsize numberNodes, int gSize, int avgInDegree, int avgNodeDegree){
	n = numberNodes;
	nlist = new Set[n];
	name = new ULI[n];
	float probability = (float(avgInDegree)) / (gSize - 1);
	Vertex p = probability * 10000;
	Vertex numberOfGroups = n / gSize;
	
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
	
	// To reach the desired overall average node degree, go through each node again and connect edges randomly to outside of group
	int half = (avgNodeDegree - avgInDegree) / 2.0;
	// If remaining degree is even
	if((avgNodeDegree - avgInDegree) % 2 == 0){
		for(Vertex i = 0; i < n; i++){
			randomEdges(i, gSize, half);
		}	
	}
	// If remaining degree is odd.
	else{
		for(Vertex i = 0; i < n/2; i++){
			randomEdges(i, gSize, half);
		}
		for(Vertex i = n/2; i < n; i++){
			randomEdges(i, gSize, half+1);
		}
	}
}



// Creates a graph with triangles. Only use this one if avgNodeDegree is only 1 more than avgInDegree. LastEdgeOut is the probability 
// (out of 10000) that the node's 1 edge going out of group will form.
void Graph::createTriGraph(Vsize numberNodes, int gSize, int avgInDegree, int avgNodeDegree, int lastEdgeOut){
	
	if (avgNodeDegree - avgInDegree != 1){
		cout << "This function is only for when one edge per node goes outside of its group\n";
		return;
	}
	
	n = numberNodes;
	nlist = new Set[n];
	name = new ULI[n];
	float probability = (float(avgInDegree)) / (gSize - 1);
	Vertex p = probability * 10000;
	Vertex numberOfGroups = n / gSize;
	
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
	
	// To get 1 edge per node to go out of that node's group, only half the nodes can attempt to create an edge. 
	int half = (avgNodeDegree - avgInDegree) / 2.0;
	for(Vertex i = n/2; i < n; i++){
		if(rand() % 10000 < lastEdgeOut)
			randomEdges(i, gSize, half+1);
	}

}



// Generate a triangle-free graph OR a random graph without groups. AvgNodeDegree should be > half the group size
void Graph::createGraph(Vsize numberNodes, int gSize, int avgNodeDegree){
	n = numberNodes;
	nlist = new Set[n];
	name = new ULI[n];
	
	// Generate graph without groups
	if(gSize == 0){
		//Set up the size of all the neighbor lists
		for(Vertex i = 0; i < n; i++){
	    	name[i] = i;
	    	nlist[i].init(avgNodeDegree * 2);
		}
		// All edges added are done randomly
		for(Vertex i = 0; i < n; i++){
			randomEdges(i, gSize, avgNodeDegree/2.0);
		}
	}
	
	else {
		int numberOfGroups = n / gSize;
		long int iterations = numberOfGroups * gSize;		// iterations = n, if n/groupSize doesn't have remainder
		float qi = avgNodeDegree - (gSize / 2.0);			// Number of edges going out of group for each node
		int setSize = (qi * 3)/2;							// Used for initial size of each node's neighbor list
		int half = qi / 2;									// Used to form edges going out of every group per node
		
		//Set up the size of all the neighbor lists
		for(Vertex i = 0; i < n; i++){
	    	name[i] = i;
	    	nlist[i].init(setSize);
		}
		
		// Add edges within each group in bipartite manner. Divide group in half. Connect each node in half to every node
		// in other half
		for(int i = 0; i < numberOfGroups; i++){
			bipartite(i*gSize, gSize);
		}
		
		// If last group is smaller than the rest, add edges using different parameters
		if(iterations != n){
			bipartite(iterations, n - iterations);
		}
		
		// Add remaining edges randomly outside group
		// if remaining degree is even
		if(int(qi) % 2 == 0){
			for(Vertex i = 0; i < n; i++){
				randomEdges(i, gSize, half);
			}
		}
		
		// if remaining degree is odd
		else{
			for(Vertex i = 0; i < (n/2); i++){
				randomEdges(i, gSize, half);
			}
			for(Vertex i = (n/2); i < n; i++){
				randomEdges(i, gSize, half+1);
			}
		}
	}
}



// Create triangle-free graph (groups) with desired avg degree and in-degree for each node. 
void Graph::TriangleFreeGraph(Vsize numberNodes, int gSize, int inDeg, int avgNodeDegree) {
	
	if (inDeg > (gSize / 2.0)){
		cout << "In degree cannot be greater than half the group size\n";
		return;
	}
	
	n = numberNodes;
	nlist = new Set[n];
	name = new ULI[n];
	
	int numberOfGroups = n / gSize;
	long int iterations = numberOfGroups * gSize;			// iterations = n if group size divides n w/o a remainder
	int edgesLeft = avgNodeDegree - inDeg;					// # of edges going out of group per node
	int halfEdgesLeft = edgesLeft / 2;
	double p = inDeg / (gSize/2.0);							// probability p you form an edge within the group in bipartite manner
	p = p * 10000;										
	int setSize = avgNodeDegree * 1.5;						// setting initial size of neighborlist of each node
	
	//Set up the size of all the neighbor lists
	for(Vertex i = 0; i < n; i++){
    	name[i] = i;
    	nlist[i].init(setSize);
	}
	
	// Add edges within group in bipartite manner with probability p
	for(int i = 0; i < numberOfGroups; i++){
		bipartiteP(i*gSize, gSize, p);
	}
	// If the last group is smaller than the rest, add edges within group in bipartite manner as is.
	if(iterations != n){
		bipartite(iterations, n - iterations);
	}
	
	// Add remaining edges for each node. These edges go outside group node belongs to.
	// if remaining degree is even
	if(int(edgesLeft) % 2 == 0){
		for(Vertex i = 0; i < n; i++){
			randomEdges(i, gSize, halfEdgesLeft);
		}
	}
	// if remaining degree is odd
	else{
		for(Vertex i = 0; i < (n/2); i++){
			randomEdges(i, gSize, halfEdgesLeft);
		}
		for(Vertex i = (n/2); i < n; i++){
			randomEdges(i, gSize, halfEdgesLeft + 1);
		}
	}
}



// Similiar to TriangleFreeGraph, but only for when one edge is going out of the group per node. ie AvgNodeDegree = inDeg + 1
// Last input 'lastEdgeOut' is the probability, out of 10000, that the node's 1 edge going out the group will form.
void Graph::TriangleFreeGraph(Vsize numberNodes, int gSize, int inDeg, int avgNodeDegree, int lastEdgeOut){
	
	if (inDeg > (gSize / 2.0)){
		cout << "In degree cannot be greater than half the group size\n";
		return;
	}
	
	if (avgNodeDegree - inDeg != 1){
		cout << "This function is only for when one edge goes outside of group per node\n";
		return;
	}
	
	n = numberNodes;
	nlist = new Set[n];
	name = new ULI[n];
	
	int numberOfGroups = n / gSize;
	long int iterations = numberOfGroups * gSize;
	int edgesLeft = avgNodeDegree - inDeg;
	int halfEdgesLeft = edgesLeft / 2;
	double p = inDeg / (gSize/2.0);
	p = p * 10000;
	int setSize = avgNodeDegree * 1.5;
	
	
	//Set up the size of all the neighbor lists
	for(Vertex i = 0; i < n; i++){
    	name[i] = i;
    	nlist[i].init(setSize);
	}
	
	for(int i = 0; i < numberOfGroups; i++){
		bipartiteP(i*gSize, gSize, p);
	}
	if(iterations != n){
		bipartite(iterations, n - iterations);
	}
	
	// Remaining degree is 1. Only half the nodes will attempt to form an edge going out of their group.
	for(Vertex i = (n/2); i < n; i++){
		if(rand() % 10000 < lastEdgeOut)
			randomEdges(i, gSize, halfEdgesLeft + 1);
	}
	
}


// Create Newman style graph
void Graph::Newman(Vertex numNodes, int avgDeg, double cluster, int groupSize) {
	
	long int numEdges = numNodes * avgDeg / 2;		// Number of edges needed to form graph with avgDeg
	long int actualNumEdges = 0;					// Counter for the number of edges formed so far
	double p;										// Probability p an edge will form between two nodes in the same group
	//from Newman paper
	double const1 = avgDeg/(groupSize - 1);
	double const2 = cluster*(groupSize - 1)/(groupSize - 2);
	double num = pow(cluster, 2) + (4*const1*const2);
	p = 0.5*(cluster + pow(num,0.5));
	int prob = p * 10000;							// Use this as the probability throughout function
	cout << p;
	
	//Checks if given parameters produce a p valude <= 1. If not, return from function.
	if(p > 1){
		cout << "p = " << p << " is greater than 1. Invalid parameters\n";
		return;
	}
	
	// Initialize graph
	n = numNodes;
	nlist = new Set[numNodes];
	name = new ULI[numNodes];
	for(Vertex i = 0; i < n; i++){
		name[i] = i;
		nlist[i].init(avgDeg * 2);
	}
	
	// Repeatedly create different groups with randomly chosen nodes until enough edges have been added to the graph.
	Vertex group[groupSize];
	Vertex node;
	while(actualNumEdges < numEdges){
 		//Choose n random nodes to fill the array
		for(int i = 0; i < groupSize; i++){
			node = rand() % numNodes;
			//Make sure every node is only in a group once. If a node is repeated, select another node.
			while(searchArray(node, group, i)){
				node = rand() % numNodes;
			}
			group[i] = node;
		}
		
		// Attempt to join an edge with probability p between each node within group once.
		for(int i = 0; i < groupSize - 1; i++){
			for(int j = i+1; j < groupSize; j++){
				if(prob < rand() % 10000){
					nlist[group[i]].Dinsert(group[j]);
					nlist[group[j]].Dinsert(group[i]);
					actualNumEdges++;
				}
			}
		}
	}
}


// Creates a triangle-free Newman graph by splitting group in half and adding edges in bipartite manner.
void Graph::BipNewman(Vertex numNodes, int avgDeg, double cluster, int groupSize) {
	
	long int numEdges = numNodes * avgDeg / 2;
	long int actualNumEdges = 0;
	double p;
	//from Newman paper
	double const1 = avgDeg/(groupSize - 1);
	double const2 = cluster*(groupSize - 1)/(groupSize - 2);
	double num = pow(cluster, 2) + (4*const1*const2);
	p = 0.5*(cluster + pow(num,0.5));
	// Find probability given that each node can only connect with nodes in other half of group. Almost double.
	double Pprime = (p * 2 * (groupSize - 1)) / groupSize;	
	int prob = Pprime * 10000;		// Use this as probability throughout function
	cout << p << " " << Pprime;
	
	if (Pprime > 1) {
		cout << "p = " << Pprime << " is greater than 1. Inavalid parameters\n";
		return;
	}
	
	// Initialize graph
	n = numNodes;
	nlist = new Set[numNodes];
	name = new ULI[numNodes];
	for(Vertex i = 0; i < n; i++){
		name[i] = i;
		nlist[i].init(avgDeg * 2);
	}
	
	// Repeatedly create different groups with randomly chosen nodes until enough edges have been added to the graph.
	Vertex group[groupSize];
	Vertex node;
	while(actualNumEdges < numEdges){
 		//Choose n random nodes to fill the array
		for(int i = 0; i < groupSize; i++){
			node = rand() % numNodes;
			//Make sure every node is only in a group once. If a node is repeated, select another node.
			while(searchArray(node, group, i)){
				node = rand() % numNodes;
			}
			group[i] = node;
		}
		
		// Attempt to create an edge with probability p between each node in one half of group with every node in second half
		for(int i = 0; i < groupSize/2; i++){
			for(int j = groupSize/2; j < groupSize; j++){
				if(Pprime < rand() % 10000){
					nlist[group[i]].Dinsert(group[j]);
					nlist[group[j]].Dinsert(group[i]);
					actualNumEdges++;
				}
			}
		}
	}
}

bool Graph::searchArray(Vertex k, Vertex array[], int size){
	for(int i = 0; i < size; i++){
		if(k == array[i])
			return true;
	}
	return false;
}

// Add edges within a group in bipartite manner. Divide group in half. Connect each node in one half
// to every node in other half.
void Graph::bipartite(Vertex startNode, int groupSize){
	Vertex midNode = startNode + (groupSize/2);
	Vertex lastNode = startNode + groupSize;
	for(Vertex i = startNode; i < midNode; i++){
		for(Vertex k = midNode; k < lastNode; k++){
			nlist[i].Dinsert(k);
			//Connect node in second half to node from first half.
			nlist[k].Dinsert(i);
		}
	}
}

// Add edges within in a group in bipartite manner. Divide group in half probability p, attempt to connect
// each node in one half with every node in other half.
void Graph::bipartiteP(Vertex startNode, int groupSize, int p){
	Vertex midNode = startNode + (groupSize/2);
	Vertex lastNode = startNode + groupSize;
	for(Vertex i = startNode; i < midNode; i++){
		for(Vertex k = midNode; k < lastNode; k++){
			if(rand() % 10000 < p){
				nlist[i].Dinsert(k);
				//Connect node in second half to node from first half.
				nlist[k].Dinsert(i);
			}
		}
	}
}

// Randomly add edges between node i and nodes outside the group node i belongs to.
void Graph::randomEdges(Vertex i, int groupSize, int half){
	Vertex randomNode;     			//The node to add.
	bool prevConnected;
	bool sameGroup;
	bool itself;
	// There are groups
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



void Graph::gnr(Vsize numberNodes, int avgDegree){
	n = numberNodes;
	nlist = new Set[n];
	name = new ULI[n];
	double p;
	p = (double(avgDegree)/(n - 1));
	p = p * 10000;

	// initialize size for each node in the graph
	for(Vertex i = 0; i < n; i++){
	    name[i] = i;
	    nlist[i].init(avgDegree * 2);
	}
	
	for(Vertex i = 0; i < n-1; i++){
		for(Vertex j = i + 1; j < n; j++){
			if(rand()%10000 < p){
				nlist[i].Dinsert(j);
				nlist[j].Dinsert(i);
			}
		}
	}
}
#endif	/* _GRAPH_HPP */