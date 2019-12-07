//Please see the header file and readme for documentation.


#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <limits.h>
#include <boost/array.hpp> 
#include "RootedGraph.h"





using namespace std;


#include "nauty26r12/nausparse.h"


vertex::vertex(int i, int color1):in(false),isIndex(false),index(i),curIndex(-1),distance(INT_MAX),color(color1),component(INT_MAX),neighbors({}){};


//compute the rooted graph of radius r centered at v, and the distances
rootedGraph::rootedGraph(vertex* v, int r1){
	r=r1;
	vertices={{v}};

	//set local data
	v->in=true; //used only in this computation 
	v->distance=0; //used in other computations
	
	vector<vertex*> curStack={v};
	int depth=0;

	//breadth-first search
	while (depth<r){
		vector<vertex*> nextStack;
		while(curStack.size()>0){
			vertex* curV=curStack.back();
			curStack.pop_back();
			for (int i=0;i<curV->neighbors.size();i++){
				if(curV->neighbors[i]->in==false){ //not seen previously
					curV->neighbors[i]->distance=depth+1;
					curV->neighbors[i]->in=true;
					nextStack.push_back(curV->neighbors[i]);
				}
			}
		}
		depth++;
		vertices.push_back(nextStack);
		curStack=nextStack;
	}

	//reset data
	for (int i=0;i<vertices.size();i++){for (int j=0;j<vertices[i].size();j++){vertices[i][j]->in=false;}}
}


//destructor: clears all local data
rootedGraph::~rootedGraph(){
	for (int i=0;i<vertices.size();i++){for (int j=0;j<vertices[i].size();j++){
		vertices[i][j]->curIndex=-1;
		vertices[i][j]->distance=INT_MAX;
		vertices[i][j]->component=INT_MAX;
	}}
	vertices={};
}

// Checks if atoms in the rooted graph satisfy the (repeated) pattern. For example, if pattern={4,2} this will return true if the atoms in shells 0, 2, 4, .. have four neighbors and atoms in shells  1,3,5,... have two neighbors.
bool rootedGraph::checkValences(vector<int> pattern){
	for (int i=0;i<vertices.size();i++){
		int desiredValence=pattern[i%pattern.size()];
		for (int j=0;j<vertices[i].size();j++){
			if (vertices[i][j]->neighbors.size()!=desiredValence){return false;}
		}
	}
	return true;
}





eClass::eClass(int type1, int r1, vector<vector<int> > data1):type(type1),r(r1),data(data1),counts({}),examples({{}}),ranks({}),freqs({}){
	//compute the hash key using the hash-combine method in boost
	key=0;
	if (type==1){
		key=(data.size())*(data.size()+1)/2;
		for (int i=0;i<data.size();i++){for (int j=i;j<data[i].size();j++){
			key=key^(data[i][j]+ 0x9e3779b9 + (key << 6) + (key >> 2));
		}}
	}
	else{
		for (int i=0;i<data.size();i++){key=key+data[i].size();}
		for (int i=0;i<data.size();i++){for (int j=0;j<data[i].size();j++){
			key=key^(data[i][j]+ 0x9e3779b9 + (key << 6) + (key >> 2));
		}}
	}
	
}
	
void eClass::resize(int numPreps){

	while (counts.size()<numPreps){counts.push_back(0);}
	while(freqs.size()<numPreps){freqs.push_back(0);}
	while (examples.size()<numPreps){examples.push_back({});}
}

bool operator==(const eClass& lhs, const eClass& rhs)
{
	if (lhs.type!=rhs.type){return false;}
	if (lhs.r!=rhs.r){return false;}
	if (lhs.key!=rhs.key){return false;}

	if (lhs.data.size()!=rhs.data.size()){return false;}

	for (int i=0;i<lhs.data.size();i++){
		if (lhs.data[i].size()!=rhs.data[i].size()){return false;}
		for (int j=0;j<lhs.data[i].size();j++){
			if (lhs.data[i][j]!=rhs.data[i][j]){return false;}

		}
	}
	return true;

}	


void eClass::print(){
	if (type==0){//see nauty documentation for format
		for (int i=0;i<4;i++){
			if (i==0){cout<<"d: ";}
			else if (i==1){cout<<"v: ";}
			else if (i==2){cout<<"e: ";}
			else if (i==3){cout<<"ptn: ";}
			for (int j=0;j<data[i].size();j++){
				cout<<data[i][j];
				if (j<data[i].size()-1){
					cout<<",";
				}
			}
			cout<<endl;
		}
			
	}
	else if (type==1){
		for (int i=0;i<data.size();i++){for (int j=i;j<data.size();j++){
			if (data[i][j]==1){
				cout<<"("<<i<<","<<j<<")";
			}
			else if (data[i][j]>1){
				cout<<data[i][j]<<"x("<<i<<","<<j<<")";
			}
			if ((j<data.size()-1) or (i<data.size()-2)){if (data[i][j]>0){
				cout<<",";
			}}
		}}
	}
	else if (type==2){

		for (int i=0;i<data[0].size();i++){
			if (data[0][i]==1){
				cout<<data[0][i]<<" "<<i+1<<"-ring";
			}
			else if (data[0][i]>1){
				cout<<data[0][i]<<" "<<i+1<<"-rings";
			}
			if ((data[0][i]>0) and (i<data[0].size()-1)){cout<<", ";}
		}
	}


	else{
		cout<<"(";
		for (int i=0;i<data.size();i++){
			cout<<data[i].size();
			bool standard=true;
			for (int j=0;j<data[i].size();j++){
				bool test=((i%2==0) and (data[i][j]!=4)) or ((i%2!=0) and (data[i][j]!=2));
				if (test){
					if (standard==true){
						cout<<"(";
						standard=false;
					}
					cout<<data[i][j]<<",";
				}
			}
			if (!standard){cout<<")";}
					
			if (i<data.size()-1){cout<<",";}
		}
		cout<<")";
	}

	

		
}

 

void eClass::print(ostream &fs){
	if (type==0){//see nauty documentation for format
		for (int i=0;i<4;i++){
			if (i==0){fs<<"d: ";}
			else if (i==1){fs<<"v: ";}
			else if (i==2){fs<<"e: ";}
			else if (i==3){fs<<"ptn: ";}
			for (int j=0;j<data[i].size();j++){
				fs<<data[i][j];
				if (j<data[i].size()-1){
					fs<<",";
				}
			}
			fs<<endl;
		}
			
	}
	else if (type==1){
		for (int i=0;i<data.size();i++){for (int j=i;j<data.size();j++){
			if (data[i][j]==1){
				fs<<"("<<i<<","<<j<<")";
			}
			else if (data[i][j]>1){
				fs<<data[i][j]<<"x("<<i<<","<<j<<")";
			}
			if ((j<data.size()-1) or (i<data.size()-2)){if (data[i][j]>0){
				fs<<",";
			}}
		}}
	}
	else if (type==2){

		for (int i=0;i<data[0].size();i++){
			if (data[0][i]==1){
				fs<<data[0][i]<<" "<<i+1<<"-ring";
			}
			else if (data[0][i]>1){
				fs<<data[0][i]<<" "<<i+1<<"-rings";
			}
			if ((data[0][i]>0) and (i<data[0].size()-1)){fs<<", ";}
		}
	}


	else{
		fs<<"(";
		for (int i=0;i<data.size();i++){
			fs<<data[i].size();
			bool standard=true;
			for (int j=0;j<data[i].size();j++){
				bool test=((i%2==0) and (data[i][j]!=4)) or ((i%2!=0) and (data[i][j]!=2));
				if (test){
					if (standard==true){
						fs<<"(";
						standard=false;
					}
					fs<<data[i][j]<<",";
				}
			}
			if (!standard){fs<<")";}
					
			if (i<data.size()-1){fs<<",";}
		}
		fs<<")";
	}

	

		
}

void eClass::printWithStats(std::ofstream& fs, int numExamples){
	this->print(fs);
	fs<<endl;
	fs<<"Frequencies: ";
	for (int j=0;j<freqs.size();j++){fs<<freqs[j]<<" ";}
	fs<<endl;


	for (int j=0;j<freqs.size();j++){
		fs<<"Examples in preparation "<<j<<": ";
		for (int k=0;k<min(10,(int) examples[j].size());k++){
			fs<<examples[j][k]<<" ";
		}
		fs<<endl;
	}
}





eClass* rootedGraph::valenceProfile(){
	vector<vector<int> > valences={};
	for (int i=0;i<=r;i++){
		vector<int> curShell={};
		for (int j=0;j<vertices[i].size();j++){curShell.push_back(vertices[i][j]->neighbors.size());}
		std::sort(curShell.begin(),curShell.end());
		valences.push_back(curShell);
	}
	return new eClass(3,r,valences);
}

eClass* rootedGraph::shellCount(){
	vector<int> shellCounts={};
	for (int i=0;i<=r;i++){shellCounts.push_back(vertices[i].size());}

	return new eClass(4,r,{shellCounts});
}



//computes the rank of the first homology group of the shell annuli of the rooted graph, using the formula rank(H1)= #components-#vertices+#edges
vector<vector<int> > rootedGraph::computeH1Counts(){

	vector<vector<int> > H1Counts(r+1,vector<int>(r+1,0));


	//shell annulus between r1 and r2
	for (int r1=0;r1<=r;r1++){
		vector<vector<vertex*> > components={};
		//Each vertex at radius r1 starts in its own component. Note that the graph is NOT assumed to be bi-partite and the number of components decreases with r2 when r1 is fixed.
		for (int i=0;i<vertices[r1].size();i++){
			components.push_back(vector<vertex*>({vertices[r1][i]}));
			vertices[r1][i]->component=i;
		}
		int nC=components.size();	
		int nE=0;
		int nV=0;	

		for (int r2=r1;r2<=r;r2++){

			nV=nV+vertices[r2].size();


			for (int i=0;i<vertices[r2].size();i++){

				vertex* curVert=vertices[r2][i];

				vector<vertex*> toCheck=curVert->neighbors;
				int ind=0;

				//account for edges between vertices in shell r2 and ones in shells r2 and r2-1
				while (ind<toCheck.size()){
					vertex* otherVert=toCheck[ind];
					ind++;

					//deal with neighbors vertices at the same radius last (*)
					//if curVert has yet to be assigned to a component, then curVert->component=INT_MAX
					if ((curVert->component>components.size()) and ((otherVert->distance==r2) and (otherVert->index>curVert->index))){ 
						toCheck.push_back(otherVert);
					}
					//Count an edge if it is contained in the shell annulus and hasn't been seen before. 
					else if ( ((r2>r1) and (otherVert->distance<r2)) or ((otherVert->distance==r2) and (otherVert->index>curVert->index))){
						nE++;
						//check if the edge kills a component
						if (curVert->component!= otherVert->component){
							//otherVert ha has not been seen (note: it must also be in shell r2)
							if (otherVert->component>components.size()){
								//our earlier step (*) ensures that curVert has been seen
								otherVert->component=curVert->component;
								components[curVert->component].push_back(otherVert);
							}
							else if (curVert->component<components.size()){//vertex previously seen, merge components
								nC--;
								vector<vertex*> otherComponent=components[otherVert->component];
								for (int k=0;k<otherComponent.size();k++){
									otherComponent[k]->component=curVert->component;
									components[curVert->component].push_back(otherComponent[k]);	
								}
								otherComponent={};
							}
							else{
								curVert->component=otherVert->component;
								components[otherVert->component].push_back(curVert);	
							}
						}
					}
				}
				
			}
			H1Counts[r1][r2]=nC-nV+nE;

		}
		//reset components
		for (int i=r1;i<=r;i++){
			for (int j=0;j<vertices[i].size();j++){vertices[i][j]->component=INT_MAX;}
		}
	}
	return H1Counts;
}



//computes the H1 Barcode using Mobius inversion.
eClass* rootedGraph::H1Barcode(vector<vector<vector<vector<int> > > > mobius){

	vector<vector<int> > intervals(r+1,vector<int>(r+1,0));
	vector<vector<int> > counts=computeH1Counts();

	for (int i=0;i<r+1;i++){for (int j=i;j<r+1;j++){
		int temp=0;
		for (int m=i;m<=j;m++){for (int n=m;n<=j;n++){
			temp=temp+counts[m][n]*mobius[i][j][m][n]; //note: Mobius function computed with opposite ordering
		}}
		intervals[i][j]=temp;

	}}

	return new eClass(1,r,intervals);
}
		


/*
Breadth-first search to find all paths between the base vertex (the source) and the sink so that the distance to the source is monotonically increasing. The code assumes that the distance from the base vertex to other vertices at distances less than or equal to d(source,sink) have already been computed.
*/
vector<vector<vertex*> > vertex::findPaths(vertex* sink, bool global){
	vector<vector<vertex*> > pathStack={{sink}}; //find paths starting at the sink
	vector<vector<vertex*> > paths={};
	while (pathStack.size()>0){
		vector<vertex*> curPath=pathStack.back();
		pathStack.pop_back();
		vertex* curV=curPath.back();
		
			for (int i=0;i<curV->neighbors.size();i++){
				vertex* nextV=curV->neighbors[i];
				if (nextV==this){ //path starting at sink, ending at source
					curPath.push_back(this);
					paths.push_back(curPath);
				}
				else if (nextV->distance<curV->distance){//no backtracking
					//with the global option, primitive rings containging a certain vertex are only computed once
					if ((!global) or ((!nextV->isIndex) or (nextV->index>index))){
						vector<vertex*> newPath=curPath;
						newPath.push_back(nextV);
						pathStack.push_back(newPath);
					}
				}

		}
	}
	return paths;
}
					
	

//Computes the distance between the base vertex and another vertex. If the distance is greater than the specified limit, returns INT_MAX.
int vertex::findDistance(vertex* other, int limit){
	if (other==this){return 0;}
	this->in=true;

	int depth=0;
	vector<vertex*> curVertices={this};
	vector<vertex*> curStack={this};

	while (depth<limit){
		vector<vertex*> nextStack;
		while(curStack.size()>0){
			vertex* curV=curStack.back();
			curStack.pop_back();


			for (int i=0;i<curV->neighbors.size();i++){
				if (curV->neighbors[i]==other){//done!
					for (int j=0;j<curVertices.size();j++){curVertices[j]->in=false;} //reset data
					return depth+1;
				}
				else if (curV->neighbors[i]->in==false){
					vertex* nextV=curV->neighbors[i];
					nextV->in=true;
					nextStack.push_back(nextV);
					curVertices.push_back(nextV);
				}
			}
			
		}
		depth++;
		curStack=nextStack;
	}

	//reset data
	for (int i=0;i<curVertices.size();i++){
		curVertices[i]->in=false;
	}
	return INT_MAX; //limit reached
}



//Compute the distances of one vertex to all other vertices in a graph. Used in the Yuan and Cormack primitive rings algorithm.
vector<int> vertex::computeDistances(int r, int totalAtoms){

	vector<int> dists(totalAtoms, INT_MAX);
	this->in=true;
	this->curIndex=0;
	this->distance=0;

	vector<vertex*> curVertices={this};
	vector<vertex*> curStack={this};
	
	dists[this->index]=0;

	//breadth first search
	int depth=0;
	while (depth<r){
		vector<vertex*> nextStack;
		while(curStack.size()>0){
			vertex* curV=curStack.back();
			curStack.pop_back();
			for (int i=0;i<curV->neighbors.size();i++){
				if(curV->neighbors[i]->in==false){
					vertex* nextV=curV->neighbors[i];
					nextV->in=true;
					nextStack.push_back(nextV);
					curVertices.push_back(nextV);


					nextV->distance=depth+1;
					dists[nextV->index]=depth+1;	
				}
			}
			
		}
		depth++;
		curStack=nextStack;
	}

	//reset data
	for (int i=0;i<curVertices.size();i++){
		curVertices[i]->in=false;
		curVertices[i]->distance=INT_MAX;
	}



	return dists;
}

//Checks if a ring is primitive, using reference lists of distances to a few vertices. See Yuan and Cormack (2002). 
bool checkPrimitiveDirected(vector<vertex*> ring, vector<vector<int> > references)
{

	int sz=ring.size();
	for (int i=0;i<ring.size()-1;i++){
		for (int j=i+1;j<ring.size();j++){
			bool ok=false;
			int ringDist=min(j-i, i+sz-j);//distance along the ring
			
			//Check if non-primitivity can be ruled out using the reference distances.
			for (int k=0;k<references.size();k++){
				int d1=references[k][ring[i]->index];
				int d2=references[k][ring[i]->index];

				if (abs(d1-d2)>=ringDist){
					ok=true; //non-primitivity would violate the triangle inequality
					k=references.size();
				}
			}

			if (!ok){
				int d=ring[i]->findDistance(ring[j],ringDist); //need to perform a more costly breadth-first search to find the distance
				if (d<ringDist){
					return false;
				}
			}
			
		}
	}
	return true;
}


//Finds a list of candidate primitive rings of length <=2*r, following the Yuan and Cormack algorithm. The idea is that if the root vertex and v are both contained in a primitive ring, then v will have two neighbors that are closer to the root than it is (even ring), or one neigbhor that is the same distance from the root (odd ring).
vector<vector<vertex*> > rootedGraph::possiblePrimitive(int r, bool global){
	vector<vector<vertex*> > rings={};
	
	for (int i1=1;i1<vertices.size();i1++){for (int i2=0;i2<vertices[i1].size();i2++){if ((!global) or ((!vertices[i1][i2]->isIndex) or (vertices[i1][i2]->index>vertices[0][0]->index))){
		vertex* curV=vertices[i1][i2];
		int d=curV->distance;
		int numSame=0;
		int numShorter=0;
		vector<vertex*> sameDistanceVertices={};
		

		//checks if curV has two neighbors that are closer to the root than it is (possible even primitive ring) or at one that is the same distance (possible odd primitive ring)
		for (int j=0;j<curV->neighbors.size();j++){
			if (curV->neighbors[j]->distance==d-1){numShorter++;}
			if (curV->neighbors[j]->distance==d){
				sameDistanceVertices.push_back(curV->neighbors[j]);
				numSame++;
			}
		}
	
		//candidate even rings: concatenate paths from curV to v that have the same length
		if (numShorter>1){
			vector<vector<vertex*> > curPaths=vertices[0][0]->findPaths(curV,global);
			for (int j=0;j<curPaths.size();j++){for (int k=j+1;k<curPaths.size();k++){
				vector<vertex*> curRing=curPaths[j];
				curRing.insert(curRing.end(),curPaths[k].rbegin()+1,curPaths[k].rend()-1); //concatenate
				rings.push_back(curRing);
			}} 
		}
		//candidate odd rings
		if (numSame>0){
			vector<vector<vertex*> > curPaths=vertices[0][0]->findPaths(curV,global);
			for (int l=0;l<sameDistanceVertices.size();l++){
				vector<vector<vertex*> > otherPaths=vertices[0][0]->findPaths(sameDistanceVertices[l],global);
				for (int j=0;j<curPaths.size();j++){for (int k=0;k<otherPaths.size();k++){
					vector<vertex*> curRing=curPaths[j];
					curRing.insert(curRing.end(),otherPaths[k].rbegin(),otherPaths[k].rend()-1); //concatenate
					rings.push_back(curRing);
				}} 
			}
		}
	}}}

	return rings;
}

eClass* rootedGraph::primitiveRingProfile(std::vector<std::vector<int> > refs){
	//computes a list of candidate primitive rings
	vector<vector<vertex*> > candidateRings=possiblePrimitive(r);

	vector<int> ringProfile={};
	for (int i=0;i<candidateRings.size();i++){if (checkPrimitiveDirected(candidateRings[i],refs)){//check if a ring is primitive
		//if a primitive ring is longer than those previously detected, increase the length of the profile
		while (candidateRings[i].size()>ringProfile.size()){ringProfile.push_back(0);} 
		ringProfile[candidateRings[i].size()-1]++;
	}}

	return new eClass(2,r,{ringProfile});

}


//Computes the Mobius function of the interval poset. For historical reasons, this is the Mobius function of the poset with the opposite ordering. 
vector<vector<vector<vector<int> > > > computeMobius(int r){
	
	// for convenience
	int rad=r+1;

 	vector<vector<vector<vector<int> > > > mobius(rad, vector<vector<vector<int> > >(rad, vector<vector<int> >(rad,vector<int>(rad,0))));

	for (int i=0;i<rad;i++){for (int j=i;j<rad;j++){
		mobius[i][j][i][j]=1;
		
	}}



	for (int i=0;i<rad;i++){for (int j=i;j<rad;j++){


		for (int d=1;d<j-i+1;d++){
			for (int m1=i;m1<=i+d;m1++){
				int n1=j-(d-(m1-i));

				if (n1>=m1){
					double temp=0;

					for (int m2=i;m2<=m1;m2++){ for (int n2=n1;n2<=j;n2++){if ((m1!=m2) or (n1!=n2)){

						temp=temp-mobius[i][j][m2][n2]; 
					}}}
					mobius[i][j][m1][n1]=temp;
				}
			}
		}
	}}
	return mobius;
							
				
}


//computes canonical form for the graph isomorphism class of radius rad, using the package nauty.
eClass* rootedGraph::canonicalForm(bool primitiveCluster)
{
	vector<vector<vertex*> > verticesByColor={{}};
	for (int i=0;i<vertices.size();i++){for (int j=0;j<vertices[i].size();j++){
		while (vertices[i][j]->color >= verticesByColor.size()){verticesByColor.push_back({});}
		verticesByColor[vertices[i][j]->color].push_back(vertices[i][j]);
		vertices[i][j]->in=true;
	
	}}

	
	
	
	int ind=0;
	int numEdges=0; //note: this is the number of DIRECTED edges (so twice the number of edges)


	for (int i=0;i<verticesByColor.size();i++){
		for (int j=0;j<verticesByColor[i].size();j++){
			for (int k=0;k<verticesByColor[i][j]->neighbors.size();k++){if (verticesByColor[i][j]->neighbors[k]->in){
				numEdges++;
			}}
			verticesByColor[i][j]->curIndex=ind;
			ind++;
			
		}
	}	
	

	//initialize nauty variables
	DYNALLSTAT(int,lab,lab_sz);
	DYNALLSTAT(int,ptn,ptn_sz);
	DYNALLSTAT(int,orbits,orbits_sz);

	static DEFAULTOPTIONS_SPARSEGRAPH(options);
	statsblk stats;
	options.digraph= FALSE;
	options.getcanon = TRUE;
	options.defaultptn = FALSE;


	SG_DECL(sg);
	SG_DECL(cg);

	

	int n=ind;
	int ne=numEdges;

	int m = SETWORDSNEEDED(n);
	nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
	DYNALLOC1(int,lab,lab_sz,n,"malloc");
	DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
	DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
	SG_ALLOC(sg,n,ne,"malloc");
	sg.nv = n; 
	sg.nde = ne; 

	ind=0;
	int edgesInd=0;

	for (int l=0;l<verticesByColor.size();l++){
		for (int j=0;j<verticesByColor[l].size();j++){
			vertex* curVert=verticesByColor[l][j];
			ptn[ind]=1; 
			lab[ind]=ind;
			int curDegree=0;
			sg.v[ind]=edgesInd;
			for (int i=0;i<curVert->neighbors.size();i++){if (curVert->neighbors[i]->in){
				curDegree++;
				sg.e[edgesInd]=curVert->neighbors[i]->curIndex;
				edgesInd++;
			}}	
		
			sg.d[ind]=curDegree;
			if (curDegree==0){sg.v[ind]=0;}
			ind++;		
		}
		ptn[ind-1]=0;//end of color
	}

	
	/*
	for (int i=0;i<n;i++){
		cout<<sg.v[i]<<" ";
	}
	cout<<endl;


	for (int i=0;i<ne;i++){
		cout<<sg.e[i]<<" ";
	}
	cout<<endl;

	*/

	int *ptn2=new int[n];
	memcpy (ptn2, ptn, n *sizeof(int));
	
	sparsenauty(&sg,lab,ptn,orbits,&options,&stats,&cg);
	sortlists_sg(&cg);



	for (int i=0;i<vertices.size();i++){for (int j=0;j<vertices[i].size();j++){
		vertices[i][j]->curIndex=0;
		vertices[i][j]->in=false;
	}}


	vector<vector<int> > data={{},{},{},{}};
	for (int i=0;i<n;i++){
		data[0].push_back((int)cg.d[i]);

		data[1].push_back((int)cg.v[i]);
		data[3].push_back((int)ptn2[i]);
	}
	for (int i=0;i<ne;i++){data[2].push_back((int)cg.e[i]);}



	return new eClass(0,r,data);


	
}
