This code will produce the genomic tracks presented in figure 5 panel a. The input data is supplied in this directory. When attempting to run this script, make sure to set the working directory to location in which the input files are saved. The input files consist of merged, normalized and scaled, bedgraphs from the human noc release experiment. Specifally, time points 0, 10 hours, 15 hours, 25 hours, and two cell cycle samples. 

ERD and LRD files are suplied for annotations. These ERD and LRD files are a merged form of the original files human ERD and LRD files. When attempting top make smaller figures, Gviz will generate additional lines of annotation if the white space between two freatures of the same annotation is not small enough to display. This creates inconcsistently sized figures that adds little information to the graphic. We generated merged ERD or LRD regions into into a single regions that were seperated by 350kb or less. This allowed us to display all ERDs and LRDs over a whole chromosome without distoring the figure. 

The ERD regions are loaded a second time as a highlighted track to help visualize thier 
boundaries with respect to signal in the tracks.

All timepoint and annotations as described above are provided in this directory.

Please not that this script will save each chromosome plot to the working directory. 

This script will produce a chromosome plot for every chromsome listed in the bedgraph files.
