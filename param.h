
/*
 * param.h for FR-HIT
 * Copyright (c) 2010-2011 Beifang Niu All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _PARAM_H_
#define _PARAM_H_

#include <string>

typedef unsigned bit32_t;
typedef unsigned char bit8_t;
typedef unsigned short bit16_t;
typedef unsigned long long bit64_t;
//24bits unit
struct bit24_t { unsigned a:24; };
// seqs<4G, length <4G
typedef bit32_t ref_id_t;
typedef bit32_t ref_loc_t;

class Param {
    public:
        Param();
        ~Param();
        void SetSeedSize(int n);
        void Set4kmerParas();
        void SetRepeat();
    public:
        int ncpu;
        int chains;
        int max_dbseq_size;
        int append_dbseq_size;
        int seed_size;
        int seed_overlap;
        int min_read_size;
        int identity;
        int global;
        int align_len;
        int band;
        int max_read_size;
        int append_read_size;
        int report_repeat_hits;
        int best_kmers;      //default 20, 4-mers
        int best_nas;        //default 24, bps needed
        int max_align_hits;  //default 0
        int maxtrys;         //-t default 20, max alignment attemps
        int global_signal;   //-g default 0, 1 for global alignment control
        int global_part;     //-q default 90 when -g=1
        int lenforstep;      //-w default 1000, for 454 long reads using 2bp q-gram index step
        int outputformat;    
        int mask;	
        double evalue; // -e default 0.001, evalue cutoff;
        std::string useful_nt;    //mask lower base except "ACGT", alphabet table
        std::string nx_nt;        //mask "NXacgt"
        std::string useful_nt_nomask; //no mask
        std::string nx_nt_nomask; //only mask "NX"

        bit8_t num_mismatch[0xffffff]; //mismatch table
        bit32_t seed_bits; //mask for seed(kmer)
};

#endif

