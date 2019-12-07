/*
See the the readme for documentation.

g++ Swatches.cpp Classification.cpp RootedGraph.cpp nauty26r12/nauty.c nauty26r12/nautil.c nauty26r12/schreier.c nauty26r12/naurng.c nauty26r12/nausparse.c -Wno-write-strings -o Swatches -std=c++0x -O2 

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

//parses a comma-delimited string
vector<string> parseString(string toParse){
	vector<string> toReturn;

	std::stringstream ss(toParse);

	while(ss.good()){
		string substring;
		getline( ss, substring, ',' );
		toReturn.push_back(substring);
	}
	return toReturn;
}

int main(int argc, char** argv) {
	//f: file (SEPARATED BY COMMAS, EACH FILE A DIFFERENT DATA PREPARATION)

	vector<string> dataFiles={};
	vector<string> loadFiles={};
	int r=3;
	int type=0;
	int selection=-1;
	int Lp=-1;
	bool KL=false;
	bool shannon=false;
	string outname="";

	
	int opt;
	while ((opt = getopt(argc,argv,"f:t:r:s:p:keo:")) != EOF)
	switch(opt)
	{
		case 'f': dataFiles=parseString(optarg); break;
		case 't': type=atoi(optarg) ; break;
		case 'r': r=atoi(optarg); break;
		case 's': selection=atoi(optarg); break;
		case 'p': Lp=atoi(optarg); break;
		case 'k': KL=true; break;
		case 'e': shannon=true; break;
		case 'o': outname=optarg; break;

		case '?': fprintf(stderr, "Usage is \n -f : for names of graphs to load \n -t: for the equivalence class type \n -r: for the radius \n -s: for the selection type \n -o: for the name of the output file \n  -p: for the exponent of the Lp norm \n -k: to compute the KL divergence \n -e: to compute the Shannon entropy. \n Please see the readme for more details.");
	}


	if (dataFiles.size()==0){
		cout<<"Please provide at least one data file. See the readme for usage information."<<endl;
		return 0;
	}

	if (r<=0){
		cout<<"Please enter a valid value of the radius. See the readme for usage information."<<endl;
		return 0;
	}

	if ((type<0) or (type>5)){
		cout<<"Please enter a valid value of the type. See the readme for usage information."<<endl;
		return 0;
	}


	cout<<endl<<"Classifying environments up to ";

	if (type==1){cout<<"H1 barcode equivalence ";}
	else if (type==2){cout<<"primitive ring profile equivalence ";}
	else if (type==3){cout<<"coordination profile equivalence ";}
	else if (type==4){cout<<"shell count equivalence ";}
	else {cout<<"graph isomorphism ";} //t=0
	cout<<"at radius "<<r<<endl<<endl;
 

	empiricalDistribution* cloth=new empiricalDistribution(type,r,selection);

	cout<<"Loading data."<<endl;

	for (int i=0;i<dataFiles.size();i++){
		cout<<"Loading file "<<i<<endl;
		network* curGraph=new network(dataFiles[i]);
		cout<<"Computing the empirical distribution for file "<<i<<endl;
		cloth->computeDistribution(curGraph);
	}
	cout<<endl<<"Computation complete."<<endl;

	if (outname==""){outname=dataFiles[0];}

	if (Lp>0){
		cloth->LNorm(Lp,outname);
		cout<<"L"<<Lp<<" norm data saved to "<<outname<<"_L"<<Lp<<".txt."<<endl<<endl;
	}
	if (shannon){
		cloth->shannonEntropy(outname);
		cout<<"Shannon entropy data saved to "<<outname<<"_shannonEntropy_unrescaled.txt."<<endl<<endl;
	}

	if (KL){
		cloth->KLDivergence(outname);
		cout<<"KL Divergence data saved to "<<outname<<"_KL.txt."<<endl<<endl;
	}

	cloth->saveData_toLoad(outname);
	cout<<"Empirical distribution saved to "<<outname<<".dat."<<endl<<endl;

	cloth->saveData_toView(outname);
	cout<<"Empirical distribution data can be viewed at "<<outname<<".txt."<<endl<<endl;


	return 0;
}

