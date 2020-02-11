//
//  fasta.h
//  Created by Bruno Nevado on 10/02/2014.
// Copyright (c) 2014 Bruno Nevado. All rights reserved.
//
#include <iostream>
#include <vector>
#include <string>

class fasta {
    std::vector < std::string > matrix;
    std::vector <std::string> names;
    std::string infile;
public:
    fasta( int num_inds, int len = 0 );
    unsigned int num_lines () const {return int ( matrix.size() );}
    unsigned int num_bases () const {return int ( matrix[0].size() );}
    int is_aligned();
    void read_fasta_file ( std::string in ) ;
    void info_to_stdout(){ std::cout << "Read " << num_lines() << " sequences from file " << infile << std::endl;  };
    void free_concatenate ( fasta & al2, int verbose = 2, char miss = 'n' );
    std::string get_infile(){ return infile; };
    std::string name_at ( int ind0  ) { return names.at(ind0); };
    void append_from_vector ( std::string in, std::string name ) { matrix.push_back(in) ; names.push_back(name); }
    void write_to_file ( std::string out, int append = 0 );

};
