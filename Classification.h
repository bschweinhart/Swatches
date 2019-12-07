// Please see the readme for documentation.

#ifndef CLASSIFICATION_H
#define CLASSIFICATION_H

#include <unordered_map>
#include <vector>
#include <string>
#include "RootedGraph.h"

struct network
{
	int dataPrep; //indicates the data preparation

	std::vector<vertex*> vertices;

	void load(std::string filename); //Loads data from the format described in the readme.

	void loadRodney(std::string filename);

	network(std::string filename):vertices({}){load(filename);};
	
	network():vertices({}){};


	//Computes distances from three well-spaced vertices to the rest of the graph. Used in the primitive ring computation.	
	std::vector<std::vector<int> > computeReferences(vertex* v1);

	//avoids redundancy in primitive ring computation, stores primitive ring profile at each vertex
	void computePrimitiveRingsGlobal(int r, std::vector<int> indices,std::vector<std::vector<int> > refs);

	~network(); 
};
	


struct empiricalDistribution{
	int type;  
        // 0: Graph Isomorphism, 1: H1 Barcode, 2: Primitive Ring Profile, 3: Coordination Profile, 4: Shell Count

	int r; //radius

	int selection; 
        //i>=0: all vertices with color=i used as roots, -1: all vertices used as roots; -2: option for selecting 
        //perfectly coordinated environments in silica (assumes silicons are color 0). -3: use custom choice by using
        //the optional "indices" arguent in computeDistribution.*/

	int numPreps;
        //Allows for data from different sources to be compared. See eClass in RootedGraph.h for information on how the
        //data is stored.

	std::vector<int> numRoots;//stores the number of atomic environments for each preparation

	std::vector<std::vector<std::vector<std::vector<int> > > > mobius; //Mobius function. Used for H1 barcode. 

	std::unordered_map<int,std::vector<eClass*> > distr;
	//Dictionary for the empirical distribution. See eClass for how the key is computed and data is stored.
	//Classes with the same key are stored in a vector.
	
	//Standard initializer
	empiricalDistribution(int type1, int r1, int selection1=0):numPreps(0),type(type1),r(r1),selection(selection1),distr({}),numRoots({}){
		if (type==1){mobius=computeMobius(r);}
	}

	//void merge(empiricalDistribution* other, bool samePreparations=true); Add later.
	//Adds data from another empirical distribution. The empirical distributions must have the same radius and type
        //(but not necessarily the same selection). If samePreparations=true, uses the same numbering for data
        //preparations. If not, assumes all data preparations are new.

	//void merge(std::string filename, bool samePreparations=true); Add later.
	//Same as the previous, but loads the other empirical distribution from file.


	empiricalDistribution(std::string filename); //Initialize from data file in the format described in readme.txt.

	~empiricalDistribution(); //Deletes all equivalence classes as well.

	void computeDistribution(network* curGraph, std::vector<int> indices={});
	//Computes the empirical probability distribution. Use the "indices" option to specify a subset of indices at 
        //which to compute local environments.

	void computePrimitiveRingDistribution_faster(network* curGraph, int dataPrep=0, std::vector<int> indices={});
	//Faster method to compute primitive ring profiles: computes all primitive rings globally, then distributes to
        //each root. 



	std::vector<eClass*> convertToVector(); //Converts the dictionary to a vector.

	void load(std::string filename);//Loads from file in the format described in the readme.

	void saveData_toLoad(std::string filename); //Saves the data to filename.dat in the format described in the readme.

	void saveData_toView(std::string filename);
	//Saves the data to filename.txt in an easily interpretable format. Prints the 10 highest ranked equivalence 
        //classes for each preparation, then the n classes maximizing (frequency in preparation i - frequency in 
        //preparation j) for all i,j. 


	std::vector<std::vector<double> > LNorm(int p=1,std::string filename="");
	//Computes the lp distance between the empirical distributions for each preparation. If the filename option is 
        //used, saves the distances in "filename_L"+p+".txt".

	std::vector<std::vector<double> > KLDivergence(std::string filename="");
	//Computes the KL divergence between the empirical distributions for each preparation, under the assumption 
	//that the true distributions are absolutely continuous. If the filename option is used, saves the data in 
	//"filename_KL.txt".
       
	std::vector<double> shannonEntropy(std::string filename="");
	//Computes the (unrescaled) Shannon entropy of the empirical distributions. If the filename option is used,
        //saves the data in "filename_shannonEntropy_unresscaled.txt".

};


void saveData_toView_fromVect(std::vector<eClass*> eVect, std::string filename, int n=10, bool defaultSort=true); 
//Saves the data in an easily interpretable format. If defaultSort=false, prints data from all equivalence classes 
//without sorting. If defaultSort=true, prints the n highest ranked equivalence classes for each preparation, then the 
//n classes maximizing (frequency in preparation i - frequency in preparation j) for all i, j. 



#endif

