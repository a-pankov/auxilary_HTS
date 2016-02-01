//
//  parse_md_aligns.cpp
//  tophat_subs_dels
//
//  Created by Alex on 1/21/14.
//  Copyright (c) 2014 Alex. All rights reserved.
//

//#include "parse_md_aligns.h"


#include <iostream>
//#include <fstream>
#include <sstream>
#include <cstdio>
//#include <climits>
//#include <tr1/unordered_map>
//#include <tr1/unordered_set>
#include <algorithm>
#include <vector>
#include <numeric>
#include <api/BamReader.h>
#include <api/BamWriter.h>


#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

struct passwd *pw = getpwuid(getuid());
std::string homedir = pw->pw_dir;


using namespace BamTools;


void create_sam_index( std::string fname){
  std::cout << "\tIndexing " + fname << "...\n";
  std::string sys_cmd = "/home/apankov/bin/samtools index " + fname;
  const char * c = sys_cmd.c_str();
  system(c);
}



int main(int argc, const char * argv[])
{
  clock_t t1 = std::clock();
  time_t tm = std::time(NULL);
  std::srand(static_cast<unsigned int>( tm ));
  std::cout << std::asctime(std::localtime(&tm)) << "\n";
  
  using std::string;
  string inputFileName, fileName_good, fileName_bad, prefix;
  
  if ( argc == 2 )
  {
    inputFileName  = argv[1] ;
    //        prefix = argv[2];
    
  }
  else {
    std::cerr << "Wrong number of arguments, BAM file must follow command invocation." << std::endl;
    exit(EXIT_FAILURE);
  }
  fileName_good = prefix + "_good.bam";
  fileName_bad = prefix + "_bad.bam";
  
  //  using namespace BamTools;
  
  //    Open files for reading and writing
  
  BamReader reader;
  if (!reader.Open(inputFileName) ) {
    std::cerr << "Cant open " + inputFileName +"\n"   ;
    exit(EXIT_FAILURE);
  }
  else{
    std::cout << "Reading in BAM file:" + inputFileName + "\n";
  }
  const SamHeader header = reader.GetHeader();
  const RefVector referrences = reader.GetReferenceData();
  
  std::vector<string> chr( reader.GetReferenceCount() );
  for ( auto i = referrences.begin(); i != referrences.end(); ++i){
    chr[reader.GetReferenceID(i->RefName)] = i->RefName;
  }
  
  
  //    BamWriter writer_good, writer_bad;
  //    if(!writer_good.Open(fileName_good, header, referrences)){
  //        std::cerr << "Could not open " + fileName_good + " BAM file for writing.\n";
  //        exit(EXIT_FAILURE);
  //    }
  
  std::cout << /* "bam.Name" << '\t'<< */ "bam.IsFirstMate()" << '\t' << "!bam.IsReverseStrand()" << '\t'<< "cigar" << '\t' << "bam.Length" << '\t' <<  "chr[bam.RefID]" << '\t' << "chromInd"  << '\t' << "md" << '\t' << "bam.QueryBases" << '\t' /* << "bam.AlignedBases"<< '\t'*/ << "sub"  << '\t' << "bam.Qualities" << '\t' <<  "readInd"  << '\n' ;
  
  BamAlignment  bam;
  unsigned long long total(0), with_nm(0);
  unsigned int nm, mdInt, readPos, ind_pos;
  string strand, md, cigar;
  char sub;
  
  while ( reader.GetNextAlignment(bam) && ++total) {
    //      if(al.GetTag("NM",  nm) ){
    //        if (nm > 0){
    //          writer_good.SaveAlignment(al) && ++with_nm;
    //        }
    //      }
    if( bam.GetTag("MD",  md) ){
      cigar = "";
      
      
      std::vector< int32_t > chromInd( bam.Length );
      std::iota(chromInd.begin(), chromInd.end(), bam.Position);
      
      std::vector< int > readInd( bam.Length );
      std::iota(readInd.begin(), readInd.end(), 0);
      
      ind_pos = 0;
      
      
      std::vector<CigarOp>::const_iterator cigItr = bam.CigarData.begin();
      std::vector<CigarOp>::const_iterator cigEnd = bam.CigarData.end();
      
      for (; cigItr != cigEnd; ++cigItr) {
        cigar += std::to_string( (long long)cigItr->Length ) + cigItr->Type;
        //        cigar += cigItr->Type;
      }
      //      std::cout << '\n';
      
      //      std::cout << std::flush;
      
      for (cigItr = bam.CigarData.begin(); cigItr != cigEnd; ++cigItr) {
        
        if ( (cigItr->Type) == 'M' ) {
          ind_pos += cigItr->Length;
        } else if( (cigItr->Type) == 'I' ){
          
          std::cout << /*bam.Name << '\t'<< */ bam.IsFirstMate() << '\t' << !bam.IsReverseStrand() << '\t'<< cigar << '\t'  << bam.Length << '\t' << chr[bam.RefID] << '\t' << chromInd[ ind_pos ]-1 << '\t' << md << '_' << "INSERT" << '\t';
          for (int i=0; i < cigItr->Length; ++i) {
            std::cout << bam.QueryBases.at( ind_pos+i );
          }
/*          std::cout << '\t';
          for (int i=0; i < cigItr->Length; ++i) {
            std::cout << bam.AlignedBases.at( ind_pos+i );
          } */
          std::cout << '\t';
          for (int i=0; i < cigItr->Length; ++i) {
            std::cout << "I";
          }
          std::cout << '\t';
          for (int i=0; i < cigItr->Length; ++i) {
            std::cout << bam.Qualities.at( ind_pos+i );
          }
          std::cout << '\t' << ind_pos ;
          for (int i=1; i < cigItr->Length; ++i) {
            std::cout << ';' << ( ind_pos+i );
          }
          std::cout << '\n';
          
          
          std::transform(readInd.begin()+ind_pos, readInd.end(), readInd.begin()+ind_pos, std::bind2nd(std::plus<int>(), cigItr->Length));
          
          ind_pos += cigItr->Length;
          
          std::transform(chromInd.begin()+ind_pos, chromInd.end(), chromInd.begin()+ind_pos, std::bind2nd(std::minus<int32_t>(), cigItr->Length));
          
        } else if( (cigItr->Type) == 'N' ||  (cigItr->Type) == 'D'){
          
          std::transform(chromInd.begin()+ind_pos, chromInd.end(), chromInd.begin()+ind_pos, std::bind2nd(std::plus<int32_t>(), cigItr->Length));
          
          //          ++ind_pos;
        }
        else{
          std::cerr << "Undefined Cigar Operation (Not M,I,N,D), output incorrect\n";
        }
      }
      
      
      //      if( cigar == "33M1D18M" ) std::cout << "here"<< '\t' << chromInd[32] << '\t' << chromInd[33] << '\t' << chromInd[34] << '\t' << chromInd[35] <<'\n';
      
      readPos =0;
      std::stringstream ss;
      //      bam.GetTag("MD",  md);
      //      std::cout << md << '\t';
      ss << md;
      ss >> mdInt;
      readPos += mdInt;
      
      //      ss_iters = 0;
      //      std::cout << '\t' << bam.Length << '\t' << bam.Position << '\t' << md;
      
      //      std::cout << (bam.QueryBases == bam.AlignedBases) << '\t' << bam.Qualities << '\n';
      while( ss >> sub ){
        
        if ( sub == '^') {
          std::cout << /*bam.Name << '\t'<< */ bam.IsFirstMate() << '\t' << !bam.IsReverseStrand() << '\t'<< cigar << '\t' << bam.Length << '\t' << chr[bam.RefID] << '\t' << chromInd[readInd[readPos]-1]+1  << '\t' << md << '\t' << "DEL" << /* '\t'  << "DEL" << */ '\t';
          
          int counter(0);
          while( !std::isdigit(ss.peek()) && ++counter ){
            ss >> sub;
            std::cout << sub;
          }
          std::cout << '\t' << ( (counter==3) ? "DELETION" : "DEL")  << '\t' << readInd[readPos]  << '\n' ;
          ss >> mdInt;
          readPos += mdInt;
          continue;
        }
        
        std::cout << /* bam.Name << '\t'<< */ bam.IsFirstMate() << '\t' << !bam.IsReverseStrand() << '\t'<< cigar << '\t' << bam.Length << '\t' <<  chr[bam.RefID] << '\t' << chromInd[ readInd[readPos]]  << '\t' << md << '\t' << bam.QueryBases.at( readInd[readPos] )<< '\t' /* << bam.AlignedBases.at( readInd[readPos] )<< '\t' */ << sub  << '\t' << bam.Qualities.at( readInd[readPos] ) << '\t' <<  readInd[readPos]  << '\n' ;
        
        ++readPos;
        //        std::cout << md << '\t';
        ss >> mdInt;
        readPos += mdInt;
      }
      
    }

  }
  
  
  
  reader.Close();
  //    writer_good.Close();
  
  //    std::cout << "\nTotal count: " << total
  //    << "\nWith Mismatch: " << with_nm <<" (" << ( (double) with_nm / (double) total ) * 100.
  //    << "% of Total)" << std::flush;
  
  //    std::cout << "\nRunning Time in secs: " << (float) (std::clock() - t1)/CLOCKS_PER_SEC << '\n' << std::flush;
  
  return 0;
}
