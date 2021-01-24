// writing on a text file
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;

class writeGraph{
public:
  int groupSize;
  int numNodes;
  int averageNodeDegree;
  //Bipartite graph degree
  int bipDeg;
  int remBipDeg;
  int lastBipDeg;
  int lastRemBipDeg;
  int *group;
  
  void createGraph(int numberNodes, int groupSize, int avgNodeDegree);
  void bipartite(const char *fname);
  void joinGroups(const char *fname);
  void replaceFile(const char *fname, const char *tempGraph);
};

void writeGraph::createGraph(int numberNodes, int gSize, int avgNodeDegree){
  numNodes = numberNodes;
  groupSize = gSize;
  averageNodeDegree = avgNodeDegree;
  bipartite("example.gph");
  joinGroups("example.gph");
  replaceFile("example.gph", "tempGraph.gph");
}

void writeGraph::bipartite(const char *fname){
  ofstream myfile(fname);
  bipDeg = groupSize/2;
  remBipDeg = groupSize - bipDeg;
  group = new int[groupSize];
  int numberOfGroups = numNodes / groupSize;  //Gets truncated if the number of nodes doesn't divide evenly by the chosen group size
  int iterations = numberOfGroups * groupSize;
  
  if(myfile.is_open()){
    myfile << numNodes << '\n';
    
    for(int i = 0; i < iterations; i++){
      //At start of new group, reset group array.
      if(i % groupSize == 0){
        int temp = i;
        for(int j = 0; j < groupSize; j++){
          group[j] = temp;
          temp++;
        }
      }
    
      //Within in a group, connect each node of first half to each node of second half
      if(i % groupSize < bipDeg){
        myfile << i << ' ' << remBipDeg << '\n';
        for(int j = bipDeg; j < groupSize; j++){
          myfile << '\t' <<  group[j] << " 0 " << "0 " << '\n'; 
        } 
      }
      //Within in a group, connect each node of second half to each node of first half
      else {
        myfile << i << ' ' << bipDeg << '\n';
        for(int j = 0; j < bipDeg; j++){
          myfile << '\t' <<  group[j] << " 0 " << "0 " << '\n'; 
        } 
      }
    }
    
    //If the last group is smaller, this will connect the nodes within the last group
    if(iterations != numNodes){
      int size = numNodes - iterations;
      cout << size << '\n';
      int *lastGroup = new int[size];
      lastBipDeg = size / 2;
      lastRemBipDeg = size - lastBipDeg;
      
      //only 1 node left. 
      if (size == 1){
        myfile << iterations << ' ' << 0 << '\n';
      }
      
      //More than 1 node left, then set array and connect the nodes within the group
      else {
        int temp = iterations;
        for(int i = 0; i < size; i++){
          lastGroup[i] = temp;
          temp++;
        }
        
        for(int i = 0; i < size; i++){
          if(i < lastBipDeg){
            myfile << lastGroup[i] << ' ' << lastRemBipDeg << '\n';
            for(int j = lastBipDeg; j < size; j++){
              myfile << '\t' <<  lastGroup[j] << " 0 " << "0 " << '\n'; 
            }
          }
          
          else {
            myfile << lastGroup[i] << ' ' << lastBipDeg << '\n';
            for(int j = 0; j < lastBipDeg; j++){
              myfile << '\t' <<  lastGroup[j] << " 0 " << "0 " << '\n'; 
            }
          }
        }
      }
    }
  }
}



void writeGraph::joinGroups(const char *fname){
  ofstream myfile("tempGraph.gph");
  ifstream ifp(fname);
  float qi = averageNodeDegree - (groupSize / 2.0);
  int range = qi*2;
  int *addNodes = new int[numNodes];
  bool *prevConnected = new bool [numNodes];
  int random;         //number of edges to add to each node;
  int randomNode;     //The node to add.
  srand(time(NULL));
  int n, u, degree, v, wt, temp;
  
  if(myfile.is_open() && ifp.is_open()){
    ifp >> n;
    myfile << n << '\n';
    
    //Determine the number of edges to add to each node
    for(int i = 0; i < numNodes; i++){
      random = rand() % (range + 1);
      addNodes[i] = random;
    }
    
    //Copy the file into the temporary graph.
    for(int i = 0; i < numNodes; i++){
      
      ifp >> u >> degree;
      myfile << u << " " << addNodes[i]+degree  << '\n';
      
      for(int j = 0; j < degree; j++){
        ifp >> v >> wt >> temp;
        myfile << '\t' << v << " " << wt << " " << temp << " " << '\n';
      }
      
      //If there needs to be edges added to a node, add them.
      if(addNodes[i] != 0){
        
        //Set all nodes in the graph to false.
        for(int j = 0; j < numNodes; j++){
          prevConnected[j] = false;
        }
        
        //Node cannot connect with itself.
        prevConnected[i] = true;
        for(int k = 0; k < addNodes[i]; k++){
          randomNode = rand() % numNodes;
          //If there is already a connection, generate another random node ID
          
          while(prevConnected[randomNode] == true){
            randomNode = rand() % numNodes;
          }
          
          prevConnected[randomNode] = true;
          //print new nodes onto file.
          myfile << '\t' << randomNode << " " << "0 0 " << '\n';
        }
      }
    }
  }
}


void writeGraph::replaceFile(const char *fname, const char *tempGraph){
  const char *original = fname;
  if( remove( fname ) != 0 )
    perror( "Error deleting file" );
  else{
    rename(tempGraph, original);
  }
}