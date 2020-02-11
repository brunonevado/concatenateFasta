//
//  main.cpp
//  concatenateFasta
//
//  Created by bnevado on 26/09/2014.
//  Copyright (c) 2014 Bruno Nevado. All rights reserved.
//

#include <iostream>
#include <fstream>

#include "fasta.h"
#include "args.h"


void help(){
    std::cout << "###################\n  concatenateFasta 16022018 \n###################" << std::endl;;
    std::cout << "Concatenate any number of fasta files" << std::endl;;
    std::cout << "Usage: concatenateFasta -files in.txt -outfile out.fas -missChar N [ -pinfo partitions.txt ]" << std::endl;
    std::cout << "-files should contain list of input files to concatenate." << std::endl;
    
    std::cout << "Input files should be aligned, and in fasta format." << std::endl;
    std::cout << "Uses sequence names, so order in each file is irrelevant." << std::endl;
    std::cout << "Missing sequences are replaced with -missChar." << std::endl;
    std::cout << "Set '-verbose 0' to turn off warnings" << std::endl;
    std::cout << "Set '-pinfo file' to write partition information file for RAxML" << std::endl;

    
}


int main(int argc, const char * argv[]) {
    
    sargs myargs;
    try{
        myargs = args::getargs(argc, argv, std::vector<std::string> {"files", "outfile","missChar"}, std::vector<std::string> {}, std::vector<std::string>  {}, std::string {"pinfo"}, std::string {"verbose"}); }
    catch (std::string e){
        std::cout << " Args failed: " << e << std::endl;
        help();
        exit(1);
    }

    
    std::string infile = myargs.args_string.at(0);
    std::string outfile= myargs.args_string.at(1);
    std::string missString = myargs.args_string.at(2);
    char missChar = missString.at(0);
    
    int verb = ( myargs.args_int_optional.size() > 0 ) ? myargs.args_int_optional.at(0) : 2;
    std::string partitions = ( myargs.args_string_optional.size() > 0 ) ? myargs.args_string_optional.at(0) : "";
    
    std::vector <std::string > files;
    std::string cline;
    std::vector < int > starts;
    std::vector < int > ends;
    
    std::ifstream fh;
    fh.open(infile);
    if( !fh.is_open()){
        std::cerr << "<concatenateFasta> ERROR: cannot open for reading infile " << infile << std::endl;
        exit(1);
    }
    while (getline(fh, cline)) {
        files.push_back(cline);
    }
    
    if(verb > 0){
    std::clog << "<concatenateFasta> Read " << files.size() << " file names to concatenate from file " << infile << std::endl;
    }
    
    fasta afasta(1);
    afasta.read_fasta_file(files.at(0));
    starts.push_back(1);
    ends.push_back(afasta.num_bases());
    
    if( afasta.is_aligned() != 0 ){
        std::cerr << "<concatenateFasta> ERROR: fasta file " << files.at(0)
        << " does not seem to be aligned!" << std::endl;
        exit(1);
    }
    
    for (unsigned int i = 1; i < files.size(); i++) {
        fasta newfasta(afasta.num_lines());
        newfasta.read_fasta_file(files.at(i));
        if( newfasta.is_aligned() != 0 ){
            std::cerr << "<concatenateFasta> ERROR: fasta file " << files.at(i)
            << " does not seem to be aligned!" << std::endl;
            exit(1);
        }
        starts.push_back(afasta.num_bases() + 1);

        afasta.free_concatenate(newfasta, verb, missChar);
        ends.push_back(afasta.num_bases());

    }
    
    if(verb > 0){
    std::clog << "<concatenateFasta> Finished concatenating files, final alignment contains " << afasta.num_bases()
    << " base-pairs. Now writing to file " << outfile << std::endl;
    }
    afasta.write_to_file(outfile);
    if( partitions.length() > 0 ){
        std::ofstream fho(partitions);
        if(! fho.is_open() ){
            std::cerr << "ERROR: cant open for writing partition file " << partitions << std::endl;
            exit(1);
        }
    
        for(int i = 0; i < starts.size(); i++){
            fho << "DNA, p" << i+1 << " = " << starts.at(i) << "-" << ends.at(i) << std::endl;
        }
    }
    
    return 0;
}
