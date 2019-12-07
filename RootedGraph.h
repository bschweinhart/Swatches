//Please see the readme for documentation.

#ifndef ROOTEDGRAPH_H
#define ROOTEDGRAPH_H

struct vertex
{
	int index;
	int color; // i.e., atomic type in a bond network, or dimension of a cell in a cell complex 
	std::vector<vertex*> neighbors; //adjacent vertices 


	//Local variables used in various computations.
	bool in;
	bool isIndex;
	int curIndex;
	int distance;
	int component; //used when computing the number of components of the shell annuli in the H1 barcode computation
	std::vector<int> primitiveRingProfile; //used in global computation of primitive ring profile


	//The following functions are used in the computation of the primitive ring profile, following the Yuan and Cormack primitive ring algorithm. 
	std::vector<std::vector<vertex*> > findPaths(vertex* sink, bool global=false); //finds all distance-monotonic paths between the base vertex and the sink.
	int findDistance(vertex* other, int limit); //Computes the distance between the base vertex and another vertex, up to a limit.
	std::vector<int> computeDistances(int r, int numAtoms); //Computes the distances of one vertex to all other vertices in a graph.
	
	vertex(int i, int color1=0);

};



//eClass: short for equivalence class
struct eClass{
	int type; // 0: graph isomorphism, 1: H1 Barcode, 2: Primitive Ring Profile, 3: Coordination Profile, 4: Shell Count
	int r;
	int key; // a key to be used in a hash table
	
	std::vector<std::vector<int> > data; 
	/*The essential information representing an equivalence class. The format is different for each type:
		0: A canonical representation of the adjacancy matrix. Using the terminology from the documentation of Nauty, data={d,v,e,ptn}.
		1: A (radius+1 x radius+1) matrix I, where I(i,j)=number of intervals of the form (i,j).
		2: A single vector {{c_1,c_2,c_3...}} where c_i is the number of primitive i-rings.
		3: A vector of vectors {v_0,v_1,v_2,...} where v_i contains the valences of the vertices in the i-th shell.
		4: A single vector {{c_1,c_2,c_3...}} where c_i is the number of atoms in the i-th shell.
	
	*/
	
		
	//The following vectors represent data that is used in the computation and comparison of empirical probability distributions in classification.cpp. 
        //These are written to allow comparison of data from multiple preparations. For each, v[i]=data for the i-th preparation. 
	     
	std::vector<int> counts; //number of vertices in the equivalence class
	std::vector<double> freqs; //count divided by number of roots
	std::vector<std::vector<int> > examples; //indices of vertices in the equivalence class
	std::vector<int> ranks; //Rank[i]=r if the equivalence class is the r+1-st most common equivalence class detected in the i-th graph. 




	void resize(int numPreps);//resize internal data structures

	void print(); //prints the data from an equivalence class to the screen in an easily interpretable format
	void print(std::ostream& fs); //prints the data from an equivalence class in an easily interpretable format
	void printWithStats(std::ofstream& fs, int numExamples=10); //prints the data from an equivalence with statistics, and indices of numExamples examples from each preparation

	friend bool operator==(const eClass& lhs, const eClass& rhs);

	eClass(int type1, int r1, std::vector<std::vector<int> > data1);
};


//Data structure for the rooted graph/swatch/local atomic environment of radius r. Initializing the rooted graph sets local variables in the included vertices that are used in the computations of the various equivalence classes. These must be reset by calling the destructor before proceeding to a computation with a different local environment. 
struct rootedGraph{
	int r;

	//vertices are stored in shells. vertices[0]={root}, vertices[1]=neighbors of root, vertices[2]=vertices at distance two from the root, etc
	std::vector<std::vector<vertex*> > vertices;
	
	//eClass* graphIsomorphsimClass();
	eClass* canonicalForm(bool primitiveCluster=false); //both options need to be implemented
	eClass* H1Barcode(std::vector<std::vector<std::vector<std::vector<int> > > > mobius);
	eClass* primitiveRingProfile(std::vector<std::vector<int> > references={});
	eClass* valenceProfile();
	eClass* shellCount();

	// Checks if atoms in the rooted graph satisfy the (repeated) pattern. For example, if pattern={4,2} this will return true if the atoms in shells 0, 2, 4, .. have four neighbors and atoms in shells  1,3,5,... have two neighbors.
	bool checkValences(std::vector<int> pattern);

	std::vector<std::vector<int> > computeH1Counts(); //used in the computation of the H1 Barcode


	std::vector<std::vector<vertex*> > possiblePrimitive(int rad, bool global=false); //Finds a list of possible primitive rings containing the root.


        //After computations with one local atomic environment are complete, it is important to call this destructor which resets local data at each vertex. 
	~rootedGraph(); 

	rootedGraph(vertex* v, int r1);

};



//Compares two equivalence classes by their rank in preparation dataPrep
struct rankCompare{
	int dataPrep;
	rankCompare(int dataPrep1=0):dataPrep(dataPrep1){};
	bool operator ()(const eClass* e1,const eClass* e2)
	{
		return (e1->counts[dataPrep]>e2->counts[dataPrep]);
	}

};

//Compares two equivalence classes based on their frequencies in two different preparations
struct differenceCompare {
	int dataPrep1;
	int dataPrep2;
	differenceCompare(int dataPrep1a, int dataPrep2a):dataPrep1(dataPrep1a),dataPrep2(dataPrep2a){};
	bool operator ()(const eClass* e1,const eClass* e2)
	{
		return (e1->freqs[dataPrep1]-e1->freqs[dataPrep2]>e2->freqs[dataPrep1]-e2->freqs[dataPrep2]);
	}

};


//Checks if a ring is primitive, using distances to reference vertices to speed computation. See Yuan and Cormack (2001).
bool checkPrimitiveDirected(std::vector<vertex*> ring, std::vector<std::vector<int> > references);



//Computes the Mobius function of the interval poset. For historical reasons, this is the Mobius function of the poset with the opposite ordering. 
std::vector<std::vector<std::vector<std::vector<int> > > > computeMobius(int r);
#endif
