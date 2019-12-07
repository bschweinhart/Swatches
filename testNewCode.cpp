/*

 g++ testNewCode.cpp Classification.cpp RootedGraph.cpp nauty26r12/nauty.c nauty26r12/nautil.c nauty26r12/schreier.c nauty26r12/naurng.c nauty26r12/nausparse.c -Wno-write-strings -o testNew -std=c++0x -O2

*/



#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <boost/array.hpp> 
#include <unordered_map>
#include <utility> 
#include <iomanip> 
#include "Classification.h"



using namespace std;



int main(int argc, char** argv) {

	int radius=6;
	int type=1;
	int selection=0;
	
	empiricalDistribution* cloth=new empiricalDistribution(type,radius,selection);

	
	//network* curGraph=new network();
		
	//string filename="Positions/R10p7/connectivity_R10p7_7.dat";
	//curGraph->loadRodney(filename);
	

	/*
	rootedGraph* rGraph=new rootedGraph(curGraph->vertices[351],radius);
	eClass* curClass=rGraph->H1Barcode(computeMobius(radius));
	curClass->print(cout);
	cout<<endl;
	*/
	double elapsed=0;
	
	clock_t t;
	for (int i=0;i<100;i++){
		network* curGraph=new network();
		
		string filename="../Positions/R10p7/connectivity_R10p7_"+to_string(i)+".dat";
		curGraph->loadRodney(filename);
		t=clock();
		//for (int j=0;j<curGraph->vertices.size();j++){if (curGraph->vertices[j]->color==0){
		//	rootedGraph* rGraph=new rootedGraph(curGraph->vertices[i],radius);
		//}}
		curGraph->dataPrep=0;
		cloth->computeDistribution(curGraph);
		elapsed=elapsed+((double) clock()-t);
		
		delete curGraph;
	}


	/*
	for (int i=0;i<100;i++){
		network* curGraph=new network();
		
		string filename="../Positions/R10p6/connectivity_R10p6_"+to_string(i)+".dat";
		curGraph->loadRodney(filename);
		t=clock();
		//for (int j=0;j<curGraph->vertices.size();j++){if (curGraph->vertices[j]->color==0){
		//	rootedGraph* rGraph=new rootedGraph(curGraph->vertices[i],radius);
		//}}
		curGraph->dataPrep=1;
		cloth->computeDistribution(curGraph);
		elapsed=elapsed+((double) clock()-t);
		
		delete curGraph;
	}


	for (int i=0;i<100;i++){
		network* curGraph=new network();
		
		string filename="../Positions/R10p5/connectivity_R10p5_"+to_string(i)+".dat";
		curGraph->loadRodney(filename);
		t=clock();
		//for (int j=0;j<curGraph->vertices.size();j++){if (curGraph->vertices[j]->color==0){
		//	rootedGraph* rGraph=new rootedGraph(curGraph->vertices[i],radius);
		//}}
		curGraph->dataPrep=2;
		cloth->computeDistribution(curGraph);
		elapsed=elapsed+((double) clock()-t);
		
		delete curGraph;
	}
	*/
	cout<<"TIME: "<<elapsed/((double) CLOCKS_PER_SEC)<<endl;

	cloth->saveData_toView("data123");
	cloth->shannonEntropy("entropyTest");


}
