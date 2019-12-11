Swatches: Local Structure Classification in Graphs (beta)

DEVELOPER: Benjamin Schweinhart (schweinhart.2@osu.edu)
DATE: December 7, 2019
LICENSE: MPL 1.1/GPL 2.0/LGPL 2.1 (see license.txt)

OVERVIEW:

This software provides methods to classify local structure in graphs via empirical probability distributions of local environments. Five different notions of equivalence for local environments are supported. It is based on "Topological Similarity of Random Cell Complexes and Applications" by B. Schweinhart, J. K. Mason, and R. D. MacPherson (2016) and "Statistical Topology of Bond Networks with Applications to Silica" by B. Schweinhart, D. Rodney, and J. K. Mason (2019). Please refer to the articles for definitions.

The software can be run from the command line as described below or by including the header files "RootedGraph.h" and "Classification.h". The header files contain documentation for all functions and classes. An example "Example.cpp" is also included, which compares the local structure of Voronoi diagrams built on different point samples. 

This is a "beta" version. Please email the author at schweinhart.2@osu.edu if you experience any errors, or have any questions/comments/requests for additional functionality.


CITATIONS:

If you use this code, cite both the 2016 and 2019 papers: 

@Unpublished{2019schweinhart,
  author = {Schweinhart, B. AND Rodney, D. AND Mason, J. K.},
  title  = {Statistical Topology of Bond Networks with Applications to Silica},
  note   = {arXiv:1910.05842},
  year   = {2019},
}

@Article{2016schweinhart,
  author    = {Schweinhart, B and Mason, J. K. and MacPherson, R. D.},
  title     = {Topological similarity of random cell complexes and applications},
  journal   = {Physical Review E},
  year      = {2016},
  volume    = {93},
  number    = {6},
  publisher = {APS},
}

If you use the option to classify neighborhoods up to graph isomorphism, also cite  "Practical Graph Isomorphism, II" by McKay and Piperno (for the Nauty graph isomorphism library):

@Article{2014mckay,
  author  = {McKay, B. D. and Piperno, A.},
  title   = {Practical graph isomorphism, II},
  journal = {Journal of Symbolic Computation},
  year    = {2014},
  volume  = {60},
  pages   = {94-112},
}





DEPENDENCIES:

C++.

Boost: https://www.boost.org/.

Nauty version 2.6R12: http://pallini.di.uniroma1.it/. Download the tar.gz file from the webpage and follow the instructions provided. If you install a later version of Nauty, everything should work if you change the name of the folder to "nauty26r12" (or modify the compilation tags and the include at the top of RootedGraph.cpp).

INSTALLATION:

After installing the dependencies, compile the command line program "Swatches.cpp" as follows:

g++ Swatches.cpp Classification.cpp RootedGraph.cpp nauty26r12/nauty.c nauty26r12/nautil.c nauty26r12/schreier.c nauty26r12/naurng.c nauty26r12/nausparse.c -Wno-write-strings -o Swatches -std=c++0x -O2


INPUT FORMAT:

A graph may be inputted using the following format. The first line gives the number of vertices and the data preparation. The data preparation determines whether data from different graphs will be combined into a single empirical distribution, or a separate one. Each of the following lines corresponds to the vertex with the vertex with index 0 on the second line, index 1 on the third line, and so on. These lines are of the form "color i1 i2 i3..." where the i1, i2, i3,... are the indices of neighboring vertices. It is unecessary to enter pairs of neighboring vertices more than once, as the software will make sure the network is symmetric. Here is an example for a bi-partite ring with six vertices:
6 0
0 1 5
1 0 2
0 1 3
1 2 4
0 3 5
1 4 0

The first line indicates that the graph contains six vertices and is of data preparation 0. The second line indicates that vertex 0 is of color 0 and is adjacent to vertices 1 and 5. The next line indicates that vertex 1 is of color 1 and is adjacent to vertices 0 and 2, and so on. A larger example is included in "voronoi_uniform_10K.cfg." 
Note that a higher dimensional cell complex may be loaded by setting the color of a cell as the dimension. 



OUTPUT FORMAT:

Two output formats are available. The function "saveData_toView" in Classification.h saves the data in a format that is easy to interpret by eye (with the extension ".txt"). The first line gives the type and radius of the data, as well as the number of different data preparations. This is followed by a list of equivalence classes, where the represented data is presented in a user-friendly fashion specific to the type of the equivalence class. It should be immeadiate to interpret, except for type=0 (graph isomorphism), for which it is shown in the sparse graph representation described in the documentation of Nauty. The frequencies of the equivalence class are then given, followed by a list of the indices of n examples in each of the preparations. If defaultSort=false, the classes are displayed in the original order. If defaultSort=true, the function prints the n highest ranked equivalence classes for each preparation, then the n maximizing (frequency in preparation i - frequency in preparation j) for all i, j. An example is included in the "voronoi_comparison.txt" file. 

The function "empiricalDistribution.saveData_toLoad" in Classification.h produces data that can be re-inputted using the appropriate initializer. Here is an example of the data in the header and for the first equivalence class of an example:
1 5 0 3
10000 10000 10000 


6
0 0 0 0 2 1 
0 0 0 0 0 1 
0 0 0 0 0 2 
0 0 0 0 0 0 
0 0 0 0 0 0 
0 0 0 0 0 0 
-
2 5 12 
-
350 45 
951 630 873 351 361 
651 204 890 21 38 685 955 966 119 486 8 771 

The first line gives the type*, radius, selection method**, and number of data preparations in the data set (1, 5, 0, and 3). The next line gives the number of root atoms in each of the preparations (10000 in each). This is followed by two blank lines. Regardless of the type, the data in a equivalence class is stored in a list of vectors (actually, a vector of vectors but I say "list" of vectors for clarity.) These vectors may be of different lengths, and the number of vectors may vary between different equivalence classes. The next line is the number of vectors in the data list, followed by one line containing the contents of each vector. The next line contains the single character "-" followed by a line with the count of the number of times the equivalence class was observed in the different preparations (2 5 12). This is followed by another line with the single character '-' and finally a line for each preparation containing the indices of root vertices in the equivalence class. Data for different equivalence classes is separated by a line contianing the characters "--".

A full example is included in the "voronoi_comparison.dat" file.

*0: Graph Isomorphism, 1: H1 Barcode, 2: Primitive Ring Profile, 3: Coordination Profile, 4: Shell Count. See "Statistical Topology of Bond Networks, With Applications to Silica" for definitions. The default is t=0.

**s=-1, uses all vertices. Non-negative integers indicate that only vertices of a certain color are to be used as roots. s=-2 is a special option for silica, where only perfectly coordinated environments are used (this assumes that silica atoms are colored 0). 




COMMAND LINE:


Usage: getopt -f fname1[,fname2,fname3...] -t type -r radius [-s rootSelection] [-o outputName] [-p LpExponent] [-k] [-e]

To use the command line option, make sure you have compiled "Swatches" as described in the installation section. Different options can be selecting by using the following flags.


-f: To be used with a string of filenames separated by commas. The names of input files in the format described above. Make sure the data preparation in each input file is specified, as it will determine whether different input files are combined into a single data file, or compared separately. 

-t: To be used with an integer between 0 and 4, which specifies the equivalence relation used to classify the local environments. 0: Graph Isomorphism, 1: H1 Barcode, 2: Primitive Ring Profile, 3: Coordination Profile, 4: Shell Count. See "Statistical Topology of Bond Networks, With Applications to Silica" for definitions. The default is t=0.

-r: To be used with a positive integer, the radius of the local environments to be classified. The default is r=3.

-s: To be used with an integer greater than -3, determines which vertices are used as roots for local environments. The default, s=-1, uses all vertices. Non-negative integers indicate that only vertices of a certain color are to be used as roots. s=-2 is a special option for silica, where only perfectly coordinated environments are used (this assumes that silica atoms are colored 0). 

-o: Specifies the name of the output files. The default is to use the first filename given with the -f flag. 

-p: To be used with a positive integer p. An option for computing the Lp norm between empirical distributions. This option requires that more than one data preparation is loaded. Saves the data to the file outname+"_L"+p+".txt". The default is to not compute thie Lp norm.

-k: Include this flag if you would like to compute the KL Divergences between empirical distribution. This option requires that more than one data preparation is loaded. The divergences are saved in the file outname+"_KL.txt". The default is to not compute the KL divergence. Note that certain assumptions must be made to compute the KL divergence (in particular, equivalence classes that do not occur in both preparations are ignored), which could result in unreasonable answers.

-e: The option for (unrescaled) Shannon entropies of the empirical distributions to be computed. The entropies are saved in the file outname+"_shannonEntropy_unrescaled.txt". The default is to not compute the Shannon entropy.

Regardless of the flags used, running Swatches always saves two data files: outname+".dat" in the format to load described in the output format section above, and outname+".txt" which data for several equivalence classes in a format that is easy to interpret by eye. The second file includes the 10 highest ranked equivalence classes for each preparation, then the 10 maximizing (frequency in preparation i - frequency in preparation j) for all i, j.


For example:

	./Swatches -f voronoi_uniform_10K.cfg -r 4 -t 1 -e 
Loads the graph in "voronoi_uniform_10K.cfg," computes the empirical probability distribtion of H1 barcodes at radius 4 centered at all vertices of the data, saves the data in "voronoi_uniform_10K.cfg.dat," prints the ten most common equivalence classes in "voronoi_uniform_10K.cfg.txt," and outputs the Shannon entropy to "voronoi_uniform_10K.cfg_shannonEntropy.txt."

	./Swatches -f voronoi_uniform_10K.cfg,voronoi_lattice_10K.cfg -r 3 -t 0 -o voronoi_comparison -p 1 
Loads the graphs in "voronoi_uniform_10K.cfg" and "voronoi_lattice_10K.txt" as data preparations 0 and 1 (as specified in those files), computes the empiricial distribution of graph isomorphism classes of radius three centered at all vertices of the data, saves data files "voronoi_comparison.dat" and "voronoi_comparison.txt," and computes the L1 norm between the two empirical probability distributions saves the data in "voronoi_comparison_L1.txt". 

voronoi_uniform_10K.cfg is the 1-skeleton of a periodic Voronoi diagram of 10,000 uniformly distributed points in a square. voronoi_lattice_10K.cfg is the 1-skeleton of a periodic Voronoi diagram of 10,000 perturbed lattice points, where Gaussian noise is added to each lattice point.


OVERVIEW OF DATA STRUCTURES:

The vertex, eClass, and rootedGraph classes are declared in RootedGraph.h.

vertex: a vertex of a graph. Used in both the rootedGraph and network classes.  
eClass: short for "equivalence class." Stores the data of an equivalence class, together with information about its occurance in different data preparations (the frequency, the count of number of occurences, and indices of vertices in the equivalence class). Also computes a hash key for use in the empiricalDistribution class.
RootedGraph: constructs the rooted graph of a given radius centered at a vertex. Includes functions to compute data for each of the equivalence classes.
 

The network and empiricalDistribution classes are declared in Classification.h.
Network: the global graph data structure. Used to input data.
empiricalDistribution: computes and stores a dictionary of eClasses detected in each preparation. The dictionary maps a key (computed in eClass) with a vector of classes sharing that key. 
