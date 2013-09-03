/*
 * frhit.h for FR-HIT
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

#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "reads.h"
#include "refseq.h"
#include "align.h"
#include "param.h"
#include "utilities.h"

std::string refseqFile; 
std::string outAlignFile; 
std::string queryReadsFile;
std::ifstream finRead;
std::ifstream finDB;
std::ofstream fout;

RefSeq ref;
Param param;
ReadClass readA;

// Recruitment
unsigned int DoReadAlign() {
    ReadAlign a;
    readA.CheckFile(finRead);
    a.ImportFileFormat(readA._file_format);
    a.SetFlag('a');
    a.DoBatch(ref,readA,finRead,fout);
    return a.n_aligned;
}

void RunProcess(void) {
    if (!queryReadsFile.empty()) {
        finRead.open(queryReadsFile.c_str());
        if (!finRead) {
            std::cerr << "failed to open file: " 
                      << queryReadsFile 
                      << std::endl;
            exit(1);
        }
    } else {
        std::cerr <<"missing query file(s)\n";
        exit(1);
    }
    std::cerr << "Read recruitment:\n" 
              << "Query: " 
              << queryReadsFile 
              << "  Reference: " 
              << refseqFile 
              << std::endl;
    fout.open(outAlignFile.c_str());
    if (!fout) {
        std::cerr << "failed to open file: " 
                  << outAlignFile 
                  << std::endl;
        exit(1);
    }
    unsigned int n_aligned(0);
    readA.InitialIndex();
    n_aligned = DoReadAlign();
    finRead.close();
    fout.close();
    std::cerr << "Total number of reads recruited: " 
              << n_aligned 
              << " (" 
              << setprecision(2) 
              << 100.0*n_aligned/readA._index 
              << "%)\nDone.\nFinished at " 
              << Curr_Time()
              << "Total time consumed:  " 
              << Cal_AllTime() 
              << " secs"
              << std::endl;
}

// Usage
void usage(void) {
    std::cerr<<"\nUsage:   fr-hit v0.7 [options]\n"
        <<"       -a   <string>   reads file, *.fasta format\n"
        <<"       -d   <string>   reference genome sequences file, *.fasta format\n"
        <<"       -o   <string>   output recruitments file\n"
        <<"       -e   <double>   e-value cutoff, default="<<param.evalue<<"\n"
        <<"       -u   <int>      mask out repeats as lower cased sequence to prevent spurious hits? 1: yes; 0: no; default="<<param.mask<<"\n"
        <<"       -f   <int>      format control for output file,0:FR-HIT format; 1:PSL fromat, default="<<param.outputformat<<"\n"
        <<"       -k   <int>      k-mer size (8<=k<=12), default="<<param.seed_size<<"\n"
        <<"       -p   <int>      k-mer overlap of index (1<=p<-k), using small overlap for longer reads(454, Sanger), default="<<param.seed_overlap<<"\n"
        <<"       -c   <int>      sequence identity threshold(%), default="<<param.identity<<"\n"
        <<"       -g   <int>      use global or local alignment? 1:global; 0:local (need -m), default="<<param.global_signal<<"\n"
        <<"       -w   <int>      minimal read length to use 2bp k-mer index step to 454 long reads, default="<<param.lenforstep<<"\n"
        <<"       -m   <int>      minimal alignment coverage control for the read (g=0), default="<<param.align_len<<"\n"
        <<"       -l   <int>      length of throw_away_reads, default="<<param.min_read_size<<"\n"
        <<"       -t   <int>      maximum number of failed alingment attempts, default="<<param.maxtrys<<"\n"
        <<"       -r   [0,N]      how to report alignment hits, 0:all; N:the best top N hits for one read, default="<<param.report_repeat_hits<<"\n"
        <<"       -n   <int>      do alignment for which chain? 0:both; 1:direct only; 2:complementary only. default="<<param.chains<<"\n"
        <<"       -b   <int>      band_width of alignment, default="<<param.band<<"\n"
        <<"       -T   [0,N]      number of threads, default 1; with 0, all CPUs will be used"<< std::endl
        <<"       -h   help\n\n"
        << std::endl;
    exit(1);
}

int mGetOptions(int rgc, char *rgv[]) {
    // options
    int i;
    for (i=1; i<rgc; i++) {
        if (rgv[i][0]!='-') return i;
        switch(rgv[i][1]) {
            case 'a': queryReadsFile = rgv[++i]; break;
            case 'd': refseqFile = rgv[++i]; break;
            case 'k': param.SetSeedSize(atoi(rgv[++i])); break;
            case 'p': param.seed_overlap = atoi(rgv[++i]); break;
            case 'o': outAlignFile = rgv[++i]; break;
            case 'u': param.mask=atoi(rgv[++i]); break;
            case 'g': param.global_signal = atoi(rgv[++i]); break;
            case 'f': param.outputformat = atoi(rgv[++i]); break;
            case 'l': param.min_read_size = atoi(rgv[++i]); break;
            case 't': param.maxtrys = atoi(rgv[++i]); break;
            case 'c': param.identity = atoi(rgv[++i]); break;
            case 'm': param.align_len = atoi(rgv[++i]); break;
            case 'w': param.lenforstep = atoi(rgv[++i]); break;
            case 'e': param.evalue = atof(rgv[++i]); break;
            case 'r': param.report_repeat_hits = atoi(rgv[++i]); break;
            case 'n': param.chains=atoi(rgv[++i]); break;
            case 'b': param.band=atoi(rgv[++i]); break;
            case 'T': param.ncpu=atoi(rgv[++i]); 
                      if (param.ncpu < 0) param.ncpu = 1;
                      if (param.ncpu > MAX_THREADS) param.ncpu = MAX_THREADS;
                      break;
            case 'h':usage(); 
            case '?':usage(); 
        }
    }
    // 4-mers filtering cutoff 
    param.Set4kmerParas();
    // Repeats
    param.SetRepeat();

    return i;
}

