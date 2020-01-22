The code in this directory was used to process raw bedgraphs into their final forms for analysis. The order of operations for using this code is as follows:

1) Starting from aligned BAM files, use bedtools genomecov to generate a bedgraph. Make sure the -bga feature is used to preserve regions with zero coverage. This generates the input file  for the rpmMitochondriaNormalization.r script.
2) In addition to a bedgraph, alignmenmt stats are also needed. To generate these, sort and index the BAM file using samtools. Then use the samtools idxstats function to generate a file of read coverage on a chromosome basis.
3) Finally The same alignment idxstats are needed for the blacklisted regions. Using the blacklist files provided, create a separate BAM of these regions using samtools. Sort, index, and idxstats to get the final alignment stats per chromosome.
4) The above three files are the input for the rpmMitochondriaNormalization.r script. This will RPM normalize and scale the data.
5) Once the data is normalized, use the the biningscript.sh to standardize the bin widths. This will reduce the size of the working files as well as put all the files on the same scale to make smoothing the data and combining the data more convenient.
6) To better visualize the rerepseq data, we recommend smoothing the data using a window that is two orders of magnitude greater than your binning width. Use the bedGraphSmoothing.r script to accomplish this.

Once the bedgraphs are normalized, scaled, bined, and smoothed, replicates can be merged using bedtools unionbedg (dont forget to preserve the 0-coverage regions using -empty). Using R, excel, or any spreadsheet based program, simply average the rows for the columns containing signal for your samples to produce the final merged files.
