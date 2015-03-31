//
//  main.cpp
//  splitStrandSpecific
//
//  Created by Alex Pankov on 8/05/13.
//  Copyright (c) 2013 Alex Pankov. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <climits>
//#include <tr1/unordered_map>
//#include <tr1/unordered_set>
#include <algorithm>
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
  unsigned int max_multihits = 20;
  
  if ( argc == 3 )
  {
    inputFileName  = argv[1] ;
    prefix = argv[2];
    
  }
  else if (argc == 4){
    inputFileName  = argv[1] ;
    prefix = argv[2];
    
    int j = atoi(argv[3]);
    if(j < 2) {
      std::cerr << "k input not valid -- Must be 2 or greater\n"   ;
      exit(EXIT_FAILURE);
    }
    max_multihits = (unsigned int)j;
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
  
  std::vector<bool> nonRandChr(reader.GetReferenceCount(), true);
  std::string s;
  for ( auto i = referrences.begin(); i != referrences.end(); ++i){
    s = i->RefName;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s.find_first_not_of("chr1234567890xy") != std::string::npos)
      nonRandChr[reader.GetReferenceID(i->RefName)] = false;
    //    std::cout << s << '\t' << randChr[reader.GetReferenceID(i->RefName)] << '\n';
  }
  //  exit(EXIT_SUCCESS);
  
  
  BamWriter writer_good, writer_bad;
  if(!writer_good.Open(fileName_good, header, referrences)){
    std::cerr << "Could not open " + fileName_good + " BAM file for writing.\n";
    exit(EXIT_FAILURE);
  }
  if(!writer_bad.Open(fileName_bad, header, referrences)){
    std::cerr << "Could not open " + fileName_bad + " BAM file for writing.\n";
    exit(EXIT_FAILURE);
  }
  
  
  
  unsigned long long fStrand(0), rStrand(0),  added_xs_rev(0), added_xs_for(0), edit_xs_rev(0), edit_xs_for(0), total(0), good_al(0), singleHits(0), diff_chrom(0), matesWrongStrand(0), same_chrom(0), refid_g0(0), mapped(0), mate_mapped(0), dup(0), fail_qc(0), bad_al(0), good_map_fStrand(0), good_map_rStrand(0), weird_pos(0);
  
  //  unsigned long long total(0), good_al(0), singleHits(0);
  unsigned int nh;
  BamAlignment  al;
  uint8_t pos=43, neg=45, xs;
  
  
  while ( reader.GetNextAlignmentCore(al) && ++total) {
    
    if ( al.RefID >= 0 && ++refid_g0 && al.IsMapped() && ++mapped && al.IsMateMapped() && ++mate_mapped && !al.IsFailedQC() && ++fail_qc && nonRandChr[al.RefID] && ++good_al ){
      
      !al.IsDuplicate() && ++dup;
      if ( al.RefID != al.MateRefID ) writer_bad.SaveAlignment(al) && ++diff_chrom;
      else{
        ++same_chrom;
        al.BuildCharData();
        al.GetTag("NH",  nh);
        (nh == 1) && ++singleHits;
        if (nh >= max_multihits) continue;
        
        if( al.IsReverseStrand() && !al.IsMateReverseStrand() ){
          if (al.IsFirstMate() && ++fStrand && (al.Position > al.MatePosition) ){
            if ( al.HasTag("XS")  ){
              al.GetTag("XS", xs);
              if (xs != pos) al.EditTag("XS", "A", pos) && ++edit_xs_for; /*&& writer_bad.SaveAlignment(al);*/
            }
            else al.AddTag("XS","A",pos) && ++added_xs_for;
            writer_good.SaveAlignment(al) && ++good_map_fStrand;
          }
          else if (al.IsSecondMate() && ++rStrand && (al.Position > al.MatePosition) ){
            if ( al.HasTag("XS")  ){
              al.GetTag("XS", xs);
              if (xs != neg) al.EditTag("XS", "A", neg) && ++edit_xs_rev; /*&& writer_bad.SaveAlignment(al);*/
            }
            else al.AddTag("XS","A",neg) && ++added_xs_rev;
            writer_good.SaveAlignment(al) && ++good_map_rStrand;
          }
          else writer_bad.SaveAlignment(al) && ++weird_pos;
        }
        else if ( !al.IsReverseStrand() && al.IsMateReverseStrand() ){
          if (al.IsFirstMate() && ++rStrand && (al.Position < al.MatePosition) ){
            if ( al.HasTag("XS")  ){
              al.GetTag("XS", xs);
              if (xs != neg) al.EditTag("XS", "A", neg) && ++edit_xs_rev; /*&& writer_bad.SaveAlignment(al);*/
            }
            else al.AddTag("XS","A",neg) && ++added_xs_rev;
            writer_good.SaveAlignment(al) && ++good_map_rStrand;
          }
          else if (al.IsSecondMate() && ++fStrand && (al.Position < al.MatePosition)){
            if ( al.HasTag("XS")  ){
              al.GetTag("XS", xs);
              if (xs != pos) al.EditTag("XS", "A", pos) && ++edit_xs_for; /*&& writer_bad.SaveAlignment(al);*/
            }
            else al.AddTag("XS","A",pos) && ++added_xs_for;
            writer_good.SaveAlignment(al) && ++good_map_fStrand;
          }
          else writer_bad.SaveAlignment(al) && ++weird_pos;
        }
        else
          writer_bad.SaveAlignment(al) && ++matesWrongStrand;
      }
    }
  }
  
  reader.Close();
  writer_good.Close();
  writer_bad.Close();
  
  std::cout << "\nTotal count: " << total
  << "\nRefID ge 0: " << refid_g0
  << "\nMapped: " << mapped
  << "\nMate Mapped: "<< mate_mapped
  << "\nDuplicates: " << dup
  << "\nFail QC: " << fail_qc
  << "\nGood Alignments: " << good_al
  << "\nMates on the same chromosome count: " << same_chrom
  << "\nMates on different chromosomes count: " << diff_chrom
  << "\nSingle Hits count: " << singleHits
  << "\nForward Strand count: " << fStrand
  << "\nGood Forward Strand count: " << good_map_fStrand
  << "\nAdded XS Tag Forward Strand count: " << added_xs_for
  << "\nEdited XS Tag Forward Strand count: " << edit_xs_for
  << "\nReverse Strand count: " << rStrand
  << "\nGood Reverse Strand count: " << good_map_rStrand
  << "\nAdded XS Tag Reverse Strand count: " << added_xs_rev
  << "\nEdited XS Tag Reverse Strand count: " << edit_xs_rev
  << "\nMates with weird Positions count: " << weird_pos
  << "\nMates with non-matching strands count: " << matesWrongStrand << "\n";
  
  std::cout << "\nRunning Time in secs: " << (float) (std::clock() - t1)/CLOCKS_PER_SEC << '\n' << "\nCreating indexes:\n";
  
  create_sam_index(fileName_good);
  create_sam_index(fileName_bad);
  
  std::cout << "\nRunning Time in secs: " << (float) (std::clock() - t1)/CLOCKS_PER_SEC << '\n';
  
  return 0;
}
