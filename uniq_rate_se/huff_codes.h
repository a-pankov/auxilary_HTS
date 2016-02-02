//
//  huff_codes.h
//  ref_index
//
//  Created by Alex P on 11/15/15.
//  Copyright (c) 2015 Alex Pankov. All rights reserved.
//

#ifndef mirna_aligner_huff_codes_h
#define mirna_aligner_huff_codes_h

#include <string>
#include <vector>
#include <functional>


bool seq_has_n(const std::string & seq){
    //  std::cout << seq << '\n';
    for (size_t start = 0; start != seq.length();  ++start) {
        if (seq[start] == 'N' ) return true;
    }
    
    return false;
}

/*****************************************************
 Coding:
 C	10
 G	11
 A	00
 T	01
 ******************************************************/
void make_bitvec_2(const std::string & seq, std::vector<bool> & boolVec){
    std::vector<bool> temp( seq.length()*2 + 1 );
    temp[0] = false;
    int last =(seq.length() + 1);
    for(int i = 1; i< last; ++i){
        int bVecIndex = (2*i)-1;
        if( seq[i-1] == 'C'){
            temp[ bVecIndex ] = true;
            temp[ bVecIndex+1 ]= false;
        }
        else if( seq[i-1] == 'G'){
            temp[ bVecIndex ] = true;
            temp[ bVecIndex+1 ]= true;
        }
        else if( seq[i-1] == 'A'){
            temp[ bVecIndex ] = false;
            temp[ bVecIndex+1 ]=false;
        }
        else if( seq[i-1] == 'T'){
            temp[ bVecIndex ] = false;
            temp[ bVecIndex+1 ]= true;
        }
    }
    boolVec=temp;
}


/*****************************************************
 Coding:
 C	11
 G	10
 A	011
 T	00
 N	010
 ******************************************************/
void make_bitvec_3(const std::string & seq, std::vector<bool> & boolVec){
    std::vector<bool> temp;
    temp.push_back( true ) ;
    int last =(seq.length() + 1);
    for(int i = 1; i< last; ++i){
        //    int bVecIndex = (2*i)-1;
        if( seq[i-1] == 'C'){
            temp.push_back(true);
            temp.push_back(true);
        }
        else if( seq[i-1] == 'G'){
            temp.push_back(true);
            temp.push_back(false);
        }
        else if( seq[i-1] == 'A'){
            temp.push_back(false);
            temp.push_back(true);
            temp.push_back(true);
        }
        else if( seq[i-1] == 'T'){
            temp.push_back(false);
            temp.push_back(false);
        }
        else if( seq[i-1] == 'N'){
            temp.push_back(false);
            temp.push_back(true);
            temp.push_back(false);
        }
    }
    boolVec = temp;
}


void convert_seq_bool(const std::string & seq, std::vector<bool> & boolVec){
    if(seq_has_n(seq)) make_bitvec_3(seq, boolVec);
    else make_bitvec_2(seq, boolVec);
}

std::string convert_bool_vec_2(const std::vector<bool> & bool_vec){
    std::string temp;
    auto itr = bool_vec.begin();
    while ((++itr) != bool_vec.end() ) {
        if (*itr) {
            ++itr;
            if (*itr) temp += 'G';
            else temp += 'C';
        }else{
            ++itr;
            if (*itr) temp += 'T';
            else temp += 'A';
        }
    }
    return temp;
}

std::string convert_bool_vec_3(const std::vector<bool> & bool_vec){
    std::string temp;
    auto itr = bool_vec.begin();
    while ((++itr) != bool_vec.end() ) {
        if (*itr) {
            ++itr;
            if (*itr) temp += 'C';
            else temp += 'G';
        }else{
            ++itr;
            if (*itr){
                ++itr;
                if(*itr) temp += 'A';
                else temp += 'N';
            }
            else temp += 'T';
        }
    }
    return temp;
}

std::string convert_bool_vec(const std::vector<bool> & bool_vec){
    return bool_vec[0] ? convert_bool_vec_3(bool_vec): convert_bool_vec_2(bool_vec);
}


/*****************************************************
 Bit vector for Fastq Sequence
 ******************************************************/

//struct bitFq{
//  std::vector<bool> fq_bit;
//
//  bool operator==(const bitFq &other) const{
//    return ( fq1_bit == other.fq1_bit );
//  }
//  bitFq(const  std::string & read){
//    convert_seq_bool(read, fq_bit);
//  }
//};

/*****************************************************
 Bit Fastq Hasher
 ******************************************************/

//struct bitFqHasher{
//  std::size_t operator()(const bitFq & k) const{
//    using std::size_t;
//    using std::hash;
//    
//    return ( hash< std::vector<bool> >() (k.fq1_bit) );
//  }
//};


#endif
