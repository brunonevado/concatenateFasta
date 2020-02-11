# concatenateFasta

concatenateFasta:  concatenates several fasta files.  
  
Author: B. Nevado  
  
Usage:  
concatenateFasta -files in.txt -outfile out.fas -missChar N [ -pinfo partitions.txt ]
    -files: text file with list of input fasta files to concatenate.  
    -outfile: output fasta file to write to.  
    -missChar: character to use for missing data (e.g. N for DNA seqs)  
    -pinfo: (optional) file to write partition information to, in RAxML-compliant format.  

Output: 1 fasta file with all sequences concatenated.  Optionally also partition information file.  
    
Notes:  
    Input files should be aligned, and in fasta format.  
    Uses sequence names, so order in each file is irrelevant.  
    Missing sequences are replaced with -missChar.  
    Set '-verbose 0' to turn off warnings  
      
Installation (Linux):  
git clone https://github.com/brunonevado/concatenateFasta.git  
cd concatenateFasta  
make  
./concatenateFasta  

Folders:
/bin : executable compiled for Linux, static  
/source: source code  
/examples: test/example files, run with concatenateFasta -outfile res.fas -files files.txt -verbose 1  -missChar  N -pinfo parts.txt  
