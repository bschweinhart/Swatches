/*
In this example, we compute the empirical distributions of radius three graph isomorphism classes for one-skeletons of periodic Voronoi diagrams built on two point samples. The first is a point sample of 10,000 points uniformly distributed in a square, the second is a point sample of 10,000 points obtained by adding Gaussian noise to points on a square lattice. The Voronoi diagrams were computed using CGAL. See the readme for additional documentation. 

To compile, use the following:

 g++ Example.cpp Classification.cpp RootedGraph.cpp nauty26r12/nauty.c nauty26r12/nautil.c nauty26r12/schreier.c nauty26r12/naurng.c nauty26r12/nausparse.c -Wno-write-strings -o swatchesExample -std=c++0x -O2

*/



#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>  //std::sort
#include "Classification.h"



using namespace std;



int main(int argc, char** argv) {

	int radius=3; 
	int type=0; //classify environments using graph isomorphism
	int selection=-1; //use all vertices as roots
	cout<<"Loading data."<<endl;
	//load data from files
	network* voronoiUniform=new network("voronoi_uniform_10K.cfg"); //The one-skeleton of a periodic Voronoi diagram on
        //10K uniformly distributed points in a cube.

	network* voronoiLattice=new network("voronoi_lattice_10K.cfg");//The one-skeleton of a periodic Voronoi diagram on 
        ///10K points on a cube that were sampled by adding Gaussian noise to the vertices of the square lattice.  
	
	//initialize empirical distribution
	empiricalDistribution* cloth=new empiricalDistribution(type,radius,selection);
	
	cout<<"Computing the empirical distributions."<<endl;
	//compute empirical distributions. They will be loaded as different data preparations, as specified in the .cfg files.
	cloth->computeDistribution(voronoiUniform);
	cloth->computeDistribution(voronoiLattice);

	
	//compute the Shannon entropy
	vector<double> shannonEntropies=cloth->shannonEntropy();
	
	//Rescale the shannon entropy by log(#roots)
	double entropy1=shannonEntropies[0]/log(voronoiUniform->vertices.size());
	double entropy2=shannonEntropies[1]/log(voronoiUniform->vertices.size());
	
	cout<<endl<<"The rescaled Shannon entropy for the empirical distribution of radius 3 local environments for a Voronoi triangulation on uniform points is: "<<entropy1<<endl;
	cout<<"The rescaled Shannon entropy for the empirical distribution of radius 3 local environments for a Voronoi triangulation on perturbed lattice points is: "<<entropy2<<endl;
	cout<<"The Voronoi tesselation on uniform points exhibits more disorder in its local topology."<<endl<<endl;


	//compute the L1 distance between the empirical distributions
	vector<vector<double> > dists=cloth->LNorm(1);

	cout<<"The L1 distance between the radius 3 Voronoi empirical distributions for uniform points and perturbed lattice points is "<<dists[0][1]<<endl;


	//convert the empirical distribution to a vector
	vector<eClass*> eVect=cloth->convertToVector();
	
	//sort the equivalences classes by their frequency in the uniform sample
	sort(eVect.begin(),eVect.end(),rankCompare(0));//we added the uniform sample first, so it is data preparation 0
	
	cout<<endl<<"The most common equivalence class in the uniform empirical distribution has the following sparse graph representation: "<<endl;
	eVect[0]->print();
	cout<<"Of the "<<voronoiUniform->vertices.size()<<" local environments in the uniform sample "<<eVect[0]->counts[0]<<" are in this equivalence class."<<endl;
	cout<<"Of the "<<voronoiUniform->vertices.size()<<" local environments in the lattice sample "<<eVect[0]->counts[1]<<" are in this equivalence class."<<endl;


	//sort the equivalences classes by their frequency in the lattice  sample
	sort(eVect.begin(),eVect.end(),rankCompare(1));//we added the lattice sample second, so it is data preparation 1
	
	cout<<endl<<"The most common equivalence class in the perturbed lattice empirical distribution has the following sparse graph representation: "<<endl;
	eVect[0]->print();
	cout<<"Of the "<<voronoiUniform->vertices.size()<<" local environments in the uniform sample "<<eVect[0]->counts[0]<<" are in this equivalence class."<<endl;
	cout<<"Of the "<<voronoiUniform->vertices.size()<<" local environments in the lattice sample "<<eVect[0]->counts[1]<<" are in this equivalence class."<<endl;
	
	//sort the equivalences classes by their frequency in the uniform sample minus the frequency in the lattice sample
	sort(eVect.begin(),eVect.end(),differenceCompare(0,1));//we added the lattice sample second, so it is data preparation 1

	cout<<endl<<"The most over-represented equivalence class in the uniform sampled compared to the perturbed lattice sample has the following sparse graph representation: "<<endl;
	eVect[0]->print();
	cout<<"Of the "<<voronoiUniform->vertices.size()<<" local environments in the uniform sample "<<eVect[0]->counts[0]<<" are in this equivalence class."<<endl;
	cout<<"Of the "<<voronoiUniform->vertices.size()<<" local environments in the lattice sample "<<eVect[0]->counts[1]<<" are in this equivalence class."<<endl;

	//save the data in an easily interpretable format
	cloth->saveData_toView("voronoi_comparison");
	cout<<endl<<"To view relevant data, see voronoi_comparison.txt."<<endl;

	//save the data in a format that can be reloaded
	cloth->saveData_toLoad("voronoi_comparison");
	cout<<"Data is stored to reload in voronoi_comparison.dat."<<endl;
	

}
