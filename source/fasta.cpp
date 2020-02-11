//
//  fasta.cpp
//  Created by Bruno Nevado on 10/02/2014.
// Copyright (c) 2014 Bruno Nevado. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <locale>

#include "fasta.h"

fasta::fasta ( int num_inds, int len ){
    if(len == 0){
        matrix.reserve(num_inds);
    }
    else{
        matrix.resize(num_inds);
        
        for(int i = 0; i < num_inds; i++)
            matrix.at(i).resize(len);
    }
}

void fasta::read_fasta_file( std::string infas){
    int cind = 0;
    matrix.clear();
    std::locale loc;
    std::string line;
    std::ifstream infile_fas (infas.c_str());
    if (infile_fas.is_open())
    {
        infile = infas;
        while ( ! infile_fas.eof() )
        {
            getline(infile_fas, line);
            if( line == ""){
                continue;
            }
            if( line[0] == '>' ){
                cind++;
                matrix.resize(cind);
                std::string name = line.substr(1);
                names.push_back(name);
            }
            else {
                matrix.at(cind-1).append( line );
            }
        }
    }
    else {
        std::cerr << "ERROR (read_fasta_file): Unable to open infile " << infas << "\n" ;
        exit(1);
    }
    infile_fas.close();
    for (unsigned int l = 0; l < matrix.size(); l++) {
        for (unsigned s = 0; s < matrix.at(l).length(); s++) {
            matrix.at(l).at(s) = tolower(matrix.at(l).at(s), loc);
        }
    }
}

// CHECK WHETHER ALL SEQS HAVE SAME LENGTH
int fasta::is_aligned () {
    unsigned int previous_len = 0;
    for(unsigned int i = 0; i < matrix.size() ; i++){
        unsigned int current_len = 0;
        for (unsigned int j = 0; j <  matrix.at(i).size() ; j++) {
            if ( matrix[i][j] != char (NULL)  ){
                current_len++;
            }
            
        }
        if( i == 0){
            previous_len = current_len;
            continue;
        }
        else if (current_len != previous_len){
            return (i);
        }
    }
    
    return 0;
}

void fasta::free_concatenate(fasta &fasta2, int verbose, char miss){
    
    // checks names instead of order
    // for sequences missing in fasta2 will add Ns (raises warning)
    // sequences present in fasta2 but not fasta1 are ignored (raises warning) <- changing this of 20th April. Sequences missing are added with Ns irrespective of order
    // input fasta files should be aligned before running (does not check for this)
    
    unsigned int original_len = this->num_bases();
    
    for ( unsigned int iline = 0; iline < this->num_lines(); iline++ ) {
        bool found = false;
        for ( unsigned int iline2 = 0; iline2 < fasta2.num_lines(); iline2++ ) {
            if( this->names.at(iline) == fasta2.names.at(iline2)  ){
                found = true;
                this->matrix.at(iline).append(fasta2.matrix.at(iline2));
            }
        }
        if(!found){
            std::string seqN(fasta2.num_bases(), miss);
            this->matrix.at(iline).append(seqN.c_str());
            if(verbose > 1){
            std::clog << "<concatenateFasta2> WARNING: sequence " << this->names.at(iline) << " not found in file " << fasta2.infile << ", sequence replaced with Ns" << std::endl;
            }
        }
        
    }
    for ( unsigned int iline2 = 0; iline2 < fasta2.num_lines(); iline2++ ) {
        bool found = false;
        for( unsigned int iline = 0; iline < this->num_lines(); iline++ ) {
            if( this->names.at(iline) == fasta2.names.at(iline2)  ){
                found = true;
            }

        }
        
        if (!found){
            //this->names.push_back( fasta2.names.at(iline2));
            std::string missing( original_len, miss);
            this->append_from_vector(missing.append( fasta2.matrix.at(iline2)), fasta2.names.at(iline2));
            if(verbose > 1){
            std::clog << "<concatenateFasta2> WARNING: sequence " << fasta2.names.at(iline2) << " present in file " << fasta2.infile << " but not in " << this->infile << ". Sequence was replaced with Ns" << std::endl;
            }
        }
    }
 

}


void fasta::write_to_file( std::string out, int append ){
    std::locale loc;
   std:: ofstream outputFile;
    
    if( append == 0 )
        outputFile.open(out.c_str());
    else
        outputFile.open(out.c_str(), std::ios::app );
    
    if( !outputFile.is_open() ){
        std::cerr << "ERROR (write_to_file): unable to open for output file " << out << "\n";
        exit(1);
    }
    if( names.size() == matrix.size() ){
        for( unsigned int i = 0; i < matrix.size(); i++){
            outputFile << ">" << names[i] << "\n";
            for (unsigned int site = 0; site < matrix.at(i).length(); site++) {
                outputFile << toupper(matrix[i][site],loc); //  matrix[i][site];
            }
            outputFile << std::endl;
        }
    }
    else{
        for( unsigned int i = 0; i < matrix.size(); i++){
            outputFile << ">i" << i+1 << "\n"; //
            for (unsigned int site = 0; site < matrix.at(i).length(); site++) {
                outputFile << toupper(matrix[i][site],loc); //  matrix[i][site];
            }
            outputFile << std::endl;
        }
    }
    outputFile.close();
    
}
