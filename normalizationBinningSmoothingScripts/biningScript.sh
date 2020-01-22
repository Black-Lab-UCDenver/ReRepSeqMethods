# This script will rebin bedgraphs to standard sized bins defined by the user.
# Simply change the input.bed file to your desired input bedgraph, and change out.bed to the name of your # desired output, and then change the bin size to your desired sizes. Notice the bin size is currently 
# 100bp. Once your desired parameters have been entered, copy/paste this into a terminal to run, or run  
# per standard .sh scripts. 



rm -f tmp.db && awk 'BEGIN {printf("create table T(C TEXT,S INT,E INT,V INT); BEGIN TRANSACTION;\n");}{printf("INSERT INTO T(C,S,E,V) VALUES(\"%s\",%s,%d,%s);\n",$1,$2,$3,$4);} END {printf("COMMIT; SELECT C,(S/100)*100 as G,((S/100)+1)*100,AVG(V) FROM T GROUP BY C,G;");}' input.bed | sqlite3 -separator $'\t' tmp.db > out.bed && rm -f tmp.db

