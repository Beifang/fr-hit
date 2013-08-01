
/*
 * param.cpp for FR-HIT
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

#include "param.h"
#include <iostream>
#include <bitset>

//Aa=0 Cc=1 Gg=2 Tt=3
bit8_t alphabet[256];
void initalphabet() {  
   //convert unknown char as 'A'
   for(int i=0; i<256; i++){
       alphabet[i] = 0;
   }
   alphabet['c']=alphabet['C']=1;
   alphabet['g']=alphabet['G']=2;
   alphabet['t']=alphabet['T']=3;
}
bit8_t rev_alphabet[256];
void initrevalphabet(){ 
   //convert unknown char as 'T'
   for(int i=0; i<256; i++){
       rev_alphabet[i] = 3;
   }
   rev_alphabet['c']=rev_alphabet['C']=2;
   rev_alphabet['g']=rev_alphabet['G']=1;
   rev_alphabet['t']=rev_alphabet['T']=0;
}
bit8_t nvfilter[256];
void initnvfilter(){ 
    //'NN' filtering
   for(int i=0; i<256; i++) nvfilter[i] = 1;
   nvfilter['a']=nvfilter['A']=0;
   nvfilter['c']=nvfilter['C']=0;
   nvfilter['g']=nvfilter['G']=0;
   nvfilter['t']=nvfilter['T']=0;
}
char chain_flag[2] = {'+', '-'};
char nt_code[4] = {'A', 'C', 'G', 'T'};
char revnt_code[4] = {'T', 'G', 'C', 'A'};

Param::Param() 
    : ncpu(1)
    , chains(0)
    , max_dbseq_size(0x1000000) // 16Mb
    , append_dbseq_size(0x1000000)
    , seed_size(11)
    , seed_bits((1<<(seed_size*2))-1)
    , seed_overlap(8) //overlap
    , min_read_size(20)
    , identity(75)
    , global(0)
    , align_len(30)
    , band(4)
    , max_read_size(5000)
    , append_read_size(5000)
    , report_repeat_hits(0)
    , mask(1)
    , useful_nt("ACGT")
    , nx_nt("NXacgt")
    , useful_nt_nomask("ACGTacgt")
    , nx_nt_nomask("NX")
    , best_kmers(20)
    , best_nas(24)
    , maxtrys(20)
    , max_align_hits(0)
    , lenforstep(1000)
    , global_signal(0)
    , global_part(90)
    , outputformat(0)
    , evalue(10)
{
    initalphabet();
    initrevalphabet();
    initnvfilter();
}

Param::~Param() {
    //xxx
}

void Param::SetSeedSize(int n) {
    seed_size = n;
    seed_bits = (1<<(seed_size*2)) - 1;
#ifdef DEBUG
    bitset <32> bs2((long)seed_bits);
    cout<<" " <<bs2<<endl;
    bs2.reset();
#endif
}

void Param::Set4kmerParas() {
    // 4-mers filtering cutoff 
    best_nas = (align_len*identity) / 100;
    best_kmers = align_len - (align_len - best_nas) * 4 - 3;
    if (best_kmers < (seed_size - 3)) {
        best_kmers = seed_size - 3;
    }
}

void Param::SetRepeat() {
    // Repeats
    if (mask == 0) {
        useful_nt = useful_nt_nomask;
        nx_nt = nx_nt_nomask;
    }
}

