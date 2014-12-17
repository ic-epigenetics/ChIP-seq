#This script downloads SRA files from NCBI. 
#It requires an object called 'allfiles' containing the latter part of the path to each .sra file you wish to download. 
#E.g.

echo "SRX499/SRX499116/SRR1202457/SRR1202457
SRX499/SRX499117/SRR1202458/SRR1202458
SRX499/SRX499123/SRR1202464/SRR1202464
SRX499/SRX499124/SRR1202465/SRR1202465" > allfiles

#!/bin/bash

baseDir="http://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/"

while read LINE
do
        wget $baseDir$LINE/$LINE.sra
done < allfiles
