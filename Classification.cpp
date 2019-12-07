//Please see the readme and header file for documentation.

#define _USE_MATH_DEFINES
#include <math.h>
#include <limits.h>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility> 
#include <iomanip> 
#include <boost/array.hpp>
 

#include "RootedGraph.h"
#include "Classification.h"

using namespace std;

network::~network(){
	for (int i=0;i<vertices.size();i++){delete vertices[i];}
}



//loads an input file in the format described in the readme. If a graph already contains vertices, this adds the new data as a disconnected component.
void network::load(string filename)
{
	ifstream file(filename);

	if (file.fail()){
		cout<<"WARNING: "<<filename<<" CANNOT BE OPENED."<<endl;
	}

	string line;
	int startInd=vertices.size(); //check if there are already vertices in the graph

	getline(file,line);
	stringstream linestream(line);
	int numVerts;
	linestream>>numVerts;
	linestream>>dataPrep;

	if (numVerts!=numVerts){cout<<"WARNING: "<<filename<<" NOT IN CORRECT FORMAT."<<endl;}
	if (dataPrep!=dataPrep){cout<<"WARNING: "<<filename<<" NOT IN CORRECT FORMAT."<<endl;}
	
	//add the vertices
	for (int i=0;i<numVerts;i++){
		vertices.push_back(new vertex(i+startInd));
	}

	//color the vertices, add the edges
	for (int i=0;i<numVerts;i++){
		getline(file,line);
		stringstream linestream2(line);
		linestream2>>vertices[i]->color;

		if (vertices[i]->color!=vertices[i]->color){cout<<"WARNING: "<<filename<<" NOT IN CORRECT FORMAT."<<endl;}
		
		int neighborInd;
		while (linestream2>>neighborInd){
			//vertices[neighborInd+startInd]->neighbors.push_back(vertices[i+startInd]);
			vertices[i+startInd]->neighbors.push_back(vertices[neighborInd+startInd]);
		}
	}
	file.close();
	
	//ensure that the network is symmetric
	for (int i=startInd;i<vertices.size();i++){
		vertex* curVert=vertices[i];
		for (int j=0;j<curVert->neighbors.size();j++){
			vertex* otherVert=curVert->neighbors[j];
			bool seen=false;
			for (int k=0;k<otherVert->neighbors.size();k++){
				if (otherVert->neighbors[k]==curVert){seen=true;}
			}
			if (!seen){
				otherVert->neighbors.push_back(curVert);
			}
		}
	}
}

/*
format:
first line: #vertices dataPrep
following lines: 
  color i1 i2 i3...
  where color is the color of the vertex, and i1,i2,i3... are the indices of the neighbors
*/



void network::loadRodney(string filename)
{
	ifstream file(filename);
	string line;
	int startInd=vertices.size();


	getline(file,line);
	stringstream linestream(line);
	int numSi;
	int numO;
	linestream>>numSi;
	linestream>>numO;


	for (int i=0;i<numSi;i++){
		vertices.push_back(new vertex(i+startInd,0));
	}
	for (int i=0;i<numO;i++){
		vertices.push_back(new vertex(i+numSi+startInd,1));
	}

	getline(file,line);
	getline(file,line);
	getline(file,line);

	for (int count=0;count<numSi;count++){
		getline(file,line);
		stringstream linestream2(line);
		int neighborInd;
		while (linestream2>>neighborInd){
			neighborInd=neighborInd+numSi+startInd;
			vertices[count+startInd]->neighbors.push_back(vertices[neighborInd]);
			vertices[neighborInd]->neighbors.push_back(vertices[count+startInd]);
		}
	}
	file.close();
}

//Computes distances from three well-spaced vertices to the rest of the graph. Used in the primitive ring computation.
vector<vector<int> > network::computeReferences(vertex* v1)
{
	vector<int> ref1=v1->computeDistances(vertices.size(),vertices.size());

	//find the furthest vertex from v1 in the same connected component	
	int maxDist=0;
	vertex* v2=v1;
	for (int i=0;i<vertices.size();i++){if ((ref1[i]<vertices.size()) and (ref1[i]>maxDist)){
		maxDist=ref1[i];
		v2=vertices[i];
	}}

	vector<int> ref2=v2->computeDistances(vertices.size(),vertices.size());

	//finds the vertex that maximizes the sums of the distances to v1 and v2 	
	maxDist=0;
	vertex* v3=v1;
	for (int i=0;i<vertices.size();i++){if (ref1[i]<vertices.size()){
		if (ref1[i]+ref2[i]>maxDist){
			maxDist=ref1[i]+ref2[i];
			v3=vertices[i];
		}
	}}
	vector<int> ref3=v3->computeDistances(vertices.size(),vertices.size());
	vector<vector<int> > refs={ref1,ref2,ref3};
	return refs;
}

void network::computePrimitiveRingsGlobal(int r, vector<int> indices, vector<vector<int> > refs){
	for (int i1=0;i1<indices.size();i1++){vertices[i1]->isIndex=true;}

	for (int i1=0;i1<indices.size();i1++){

		int i=indices[i1];
		rootedGraph* rGraph=new rootedGraph(vertices[i],r);
		//eClass* curClass=rGraph->primitiveRingProfile(refs);
		vector<vector<vertex*> > candidateRings=rGraph->possiblePrimitive(r,true);
		for (int j=0;j<candidateRings.size();j++){if (checkPrimitiveDirected(candidateRings[j],refs)){//check if a ring is primitive
			//add length of primitive ring to profiles of each vertex contained in it
			for (int k=0;k<candidateRings[j].size();k++){
				vertex* curV=candidateRings[j][k];
				
				//if a primitive ring is longer than those previously detected, increase the length of the profile
				while (candidateRings[j].size()>curV->primitiveRingProfile.size()){curV->primitiveRingProfile.push_back(0);} 
				curV->primitiveRingProfile[candidateRings[j].size()-1]++;
			}
		}}
	}

}

empiricalDistribution::~empiricalDistribution()
{
	for (pair<int,vector<eClass*> > elt : distr){
		for (int i=0;i<elt.second.size();i++){
			delete elt.second[i];
		}
		elt.second={};
	}
	distr={};

}	

void empiricalDistribution::computeDistribution(network* curGraph, vector<int> indices){

	int dataPrep=curGraph->dataPrep;

	//If the data type has not been seen previously, resizes data structures in each eClass
	if (dataPrep>=numPreps){
		numPreps=dataPrep+1;
		while (numRoots.size()<dataPrep+1){numRoots.push_back(0);}

		for (pair<int,vector<eClass*> > elt : distr){for (int i=0;i<elt.second.size();i++){
			elt.second[i]->resize(numPreps);
		}}
	}

	

	//Computes the indices of vertices used as root atoms, based on the selection parameter. See Classification.h.
	if (selection!=-3){for (int i=0;i<curGraph->vertices.size();i++){
		if (selection==-2){
			rootedGraph* rGraph = new rootedGraph(curGraph->vertices[i],r);
			if (rGraph->checkValences({4,2})){indices.push_back(i);}
			delete rGraph;
		}
		else if (selection>=0){if (curGraph->vertices[i]->color==selection){indices.push_back(i);}}
		else{indices.push_back(i);}
	}}

	numRoots[dataPrep]+=indices.size();

	if (indices.size()==0){cout<<"WARNING: NO ROOT VERTICES SELECTED"<<endl;}

	//Primitive ring profile: compute reference distance matrices
	vector<vector<int> > references={};
	if (type==2){
		references=curGraph->computeReferences(curGraph->vertices[0]);


		//computes globally

		curGraph->computePrimitiveRingsGlobal(r,indices,references);
		
	}



	for (int i1=0;i1<indices.size();i1++){
		int i=indices[i1];

		//compute the rooted graph
		rootedGraph* rGraph=new rootedGraph(curGraph->vertices[i],r);


		//find the equivalence class of the rooted graph
		eClass* curClass;
		if (type==0){curClass=rGraph->canonicalForm();}
		if (type==1){curClass=rGraph->H1Barcode(mobius);}
		else if (type==2){curClass = new eClass(2,r,{curGraph->vertices[i]->primitiveRingProfile});}
		//else if (type==2){curClass=rGraph->primitiveRingProfile(references);}
		else if (type==3){curClass=rGraph->valenceProfile();}
		else if (type==4){curClass=rGraph->shellCount();}

		delete rGraph;

		//determine whether the equivalence class has been previously detected
		bool found=false;
		if (distr.find(curClass->key)!=distr.end()){//key has been seen before
			vector<eClass*> curCompare=distr[curClass->key];
			for (int j=0;j<curCompare.size();j++){
				if(*curCompare[j]==*curClass){//Equivalence class previosly detected. Update the count and the example list. 
					found=true;
					curCompare[j]->counts[dataPrep]++;
					curCompare[j]->examples[dataPrep].push_back(i);
					delete curClass;
				}
			}		
		}
		if (!found){//new equivalence class
			curClass->resize(numPreps);
			curClass->counts[dataPrep]=1;
			curClass->examples[dataPrep]={i};

			if (distr.find(curClass->key)!=distr.end()){//key has been seen before
				distr[curClass->key].push_back(curClass);
			}
			else{distr[curClass->key]={curClass};} // new key
		}
	}

	//compute the frequencies

	for (pair<int,vector<eClass*> > elt : distr){for (int i=0;i<elt.second.size();i++){
		eClass* curClass=elt.second[i];
		curClass->freqs[dataPrep]=((double) curClass->counts[dataPrep])/((double) numRoots[dataPrep]);
	}}
	
}

	
//converts a dictionary to a vector 
vector<eClass*> empiricalDistribution::convertToVector(){
	vector<eClass*> eVect={};
	for(unordered_map<int,vector<eClass*> >::iterator iter=distr.begin();iter!=distr.end();iter++){
		for (int j=0;j<(*iter).second.size();j++){
			eVect.push_back((*iter).second[j]);
		}
	}
	return eVect;
}


void saveData_toView_fromVect(vector<eClass*> eVect, string filename, int n, bool defaultSort)
{
	if (eVect.size()==0){return;}
	filename=filename+".txt";
	int numPreps=eVect[0]->freqs.size();
	ofstream fs(filename);
	fs << fixed << showpoint;
	fs << setprecision(4);

	fs<<"Type = "<<eVect[0]->type<<"  Radius="<<eVect[0]->r<<"  Number of Data Preparations="<<eVect[0]->freqs.size()<<endl; 
	fs<<endl<<"----------------------------------------------------------------------"<<endl<<endl;
	if (!defaultSort){for (int i=0;i<min(n,(int) eVect.size());i++){
		fs<<"Equivalence Class "<<i<<endl;
		eVect[i]->printWithStats(fs,10);
		fs<<endl<<endl;
	}}
	else{
		for (int j=0;j<numPreps;j++){
			fs<<"SORTED BY FREQUENCY IN PREPARATION "<<j<<endl<<endl;
			sort(eVect.begin(),eVect.end(),rankCompare(j));
			for (int i=0;i<min(n,(int) eVect.size());i++){
				fs<<"Equivalence Class "<<i<<endl;
				eVect[i]->printWithStats(fs,10);
				fs<<endl<<endl;
			}
			fs<<"----------------------------------------------------------------------";
			fs<<endl;
			fs<<endl;
		}

		for (int j=0;j<numPreps;j++){for (int k=0;k<numPreps;k++){if (k!=j){
			fs<<"SORTED BY FREQUENCY IN PREPARATION "<<j<<" MINUS FREQUENCY IN PREPARATION "<<k<<endl;
			sort(eVect.begin(),eVect.end(),differenceCompare(j,k));
			for (int i=0;i<min(n,(int) eVect.size());i++){
				fs<<"Equivalence Class "<<i<<endl;
				eVect[i]->printWithStats(fs,10);
				fs<<endl<<endl;
			}
			fs<<"----------------------------------------------------------------------";
			fs<<endl;
			fs<<endl;
		}}}
	}

	fs.close();
}

void empiricalDistribution::saveData_toView(std::string filename){
	saveData_toView_fromVect(this->convertToVector(),filename,10,true);
}
void empiricalDistribution::saveData_toLoad(std::string filename)
{

	filename=filename+".dat";
	ofstream fs(filename);
	fs<<type<<" "<<r<<" "<<selection<<" "<<numPreps<<endl;
	for (int j=0;j<numRoots.size();j++){fs<<numRoots[j]<<" ";}
	fs<<endl<<endl<<endl;

	vector<eClass*> eVect=this->convertToVector();
	if (eVect.size()==0){return;}
	for (int i=0;i<eVect.size();i++){
		fs<<eVect[i]->data.size()<<endl;
		for (int j=0;j<eVect[i]->data.size();j++){
			for (int k=0;k<eVect[i]->data[j].size();k++){fs<<eVect[i]->data[j][k]<<" ";}
			fs<<endl;
		}
		fs<<"-"<<endl;
		for (int j=0;j<eVect[i]->counts.size();j++){fs<<eVect[i]->counts[j]<<" ";}
		fs<<endl<<"-"<<endl;
		for (int j=0;j<eVect[i]->examples.size();j++){
			for (int k=0;k<eVect[i]->examples[j].size();k++){
				fs<<eVect[i]->examples[j][k]<<" ";
			}
			fs<<endl;
		}
		fs<<"--"<<endl;
		
	}
	fs.close();
}


empiricalDistribution::empiricalDistribution(std::string filename)
{
	distr={};
	ifstream file(filename);
	string line;

	getline(file,line);
	stringstream linestream(line);
	linestream>>type>>r>>selection>>numPreps;
	if (type==1){mobius=computeMobius(r);}
	
	getline(file,line);
	linestream.clear();
	linestream.str(line);
	int x;
	numRoots={};
	while (linestream>>x){numRoots.push_back(x);}

	getline(file,line);
	getline(file,line);
	while (getline(file,line)){
		int dataLength;
		linestream.clear();
		linestream.str(line);
		linestream>>dataLength;
		vector<vector<int> > data={};
		for (int i=0;i<dataLength;i++){
			getline(file,line);
			linestream.clear();
		linestream.str(line);
			vector<int> curData={};
			while (linestream>>x){curData.push_back(x);}
			data.push_back(curData);
		}
		eClass* curClass=new eClass(type,r,data);
		curClass->resize(numPreps);
		getline(file,line);
		getline(file,line);
		linestream.clear();
		linestream.str(line);
		for (int j=0;j<numPreps;j++){		
			linestream>>x;
			curClass->freqs[j]=((double)x)/numRoots[j];
			curClass->counts[j]=x;
		}
		getline(file,line);
		for (int j=0;j<numPreps;j++){
			getline(file,line);
			linestream.clear();
			linestream.str(line);
			while (linestream>>x){curClass->examples[j].push_back(x);}
		}
		getline(file,line);

		if (distr.find(curClass->key)!=distr.end()){//key has been seen before
			distr[curClass->key].push_back(curClass);
		}
		else{distr[curClass->key]={curClass};} // new key
	}

	
	file.close();
}

	

vector<vector<double> > empiricalDistribution::LNorm(int p, string filename){
	if (numPreps==1){cout<<"ONLY ONE PREPARATION - CANNOT COMPUTE L NORM"<<endl;}
	filename=filename+"_L"+to_string(p)+".txt";
	vector<eClass*> eVect=this->convertToVector();
	vector<vector<double> > dists(numPreps,vector<double>(numPreps,0.0));
	for (int i =0;i<eVect.size();i++){
		eClass* curClass=eVect[i];
		for (int j=0;j<numPreps;j++){for (int k=j+1;k<numPreps;k++){	
			double curDist=pow(fabs(curClass->freqs[j]-curClass->freqs[k]),p);
			dists[j][k]+=curDist;
		}}
	}
	for (int j=0;j<numPreps;j++){for (int k=j+1;k<numPreps;k++){
		dists[j][k]=pow(dists[j][k],1.0/p);
		dists[k][j]=dists[j][k];
	}}
	if (filename!=""){
		filename=filename+".txt";
		ofstream fs(filename);
		for (int j=0;j<numPreps;j++){
			for (int k=0;k<numPreps;k++){
				fs<<dists[j][k]<<" ";
			}
			fs<<endl;
		}
		fs.close();
	}
	
	return dists;
}

vector<vector<double> > empiricalDistribution::KLDivergence(string filename){
	if (numPreps==1){cout<<"ONLY ONE PREPARATION - CANNOT COMPUTE KL DIVERGENCES"<<endl;}
	filename=filename+"_KL.txt";
	vector<eClass*> eVect=this->convertToVector();
	vector<vector<double> > div(numPreps,vector<double>(numPreps,0.0));
	for (int i=0;i<eVect.size();i++){
		eClass* curClass=eVect[i];
		for (int j=0;j<numPreps;j++){for (int k=j+1;k<numPreps;k++){
			double p=curClass->freqs[j];
			double q=curClass->freqs[k];
			if ((p!=0) and (q!=0)){//Ignore these cases to get reasonable data!
				div[j][k]+=p*log(p/q);
				div[k][j]+=q*log(q/p);

			}
		}}
	}
	if (filename!=""){
		filename=filename+".txt";
		ofstream fs(filename);
		for (int j=0;j<numPreps;j++){
			for (int k=0;k<numPreps;k++){
				fs<<div[j][k]<<" ";
			}
			fs<<endl;
		}
		fs.close();
	}
	return div;
}


vector<double> empiricalDistribution::shannonEntropy(string filename){
	vector<eClass*> eVect=this->convertToVector();
	vector<double> ent(numPreps,0.0);
	for (int i=0;i<eVect.size();i++){for (int j=0;j<numPreps;j++){
		if (eVect[i]->counts[j]!=0){
			ent[j]+= -(eVect[i]->freqs[j])*log(eVect[i]->freqs[j]);
		}
	}}
	if (filename!=""){
		filename=filename+"_shannonEntropy_unrescaled.txt";
		ofstream fs(filename);
		for (int j=0;j<numPreps;j++){fs<<ent[j]<<" ";}
		fs<<endl;
		fs.close();
	}
	return ent;

}
