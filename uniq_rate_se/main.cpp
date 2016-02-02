//
//  main.cpp
//  uniq_rate_se
//
//  Created by Alex P on 1/27/16.
//  Copyright (c) 2016 Alex Pankov. All rights reserved.
//

#include <ctime>
#include <unistd.h>
#include <sys/types.h>

#include <zlib.h>
#include "kseq.h"
#include "huff_codes.h"
#include <iostream>
#include <unordered_set>

KSEQ_INIT(gzFile, gzread)

int main(int argc, const char * argv[]) {
    clock_t t1 = std::clock();
    time_t tm = std::time(NULL);
    std::srand(static_cast<unsigned int>( tm ));
    std::cerr << std::asctime(std::localtime(&tm));
    
    unsigned long long count(0), total(0);
    
    std::string library;

  if ( argc == 2 )
  {
    library = argv[1];
  }

    std::unordered_set<std::vector<bool>> fq_seen;
    
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(library.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0 ) {
	++total;
        auto my_seq = std::string(seq->seq.s);
        std::vector<bool> fq_bit;
        convert_seq_bool(my_seq, fq_bit);
        auto fq_search = fq_seen.find(fq_bit);
        if (fq_search == fq_seen.end()) {
            fq_seen.insert(fq_bit);
            ++count;
        }
    }
    kseq_destroy(seq);
    gzclose(fp);

    std::cerr << "Finished Mapping in " << (float) (std::clock() - t1)/CLOCKS_PER_SEC  << ".\nUnique Reads: " <<  count << "\ntotal Reads: " << total << '\n';

    
    return 0;
}
