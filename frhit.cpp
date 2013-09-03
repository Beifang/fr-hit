/*
 * frhit.cpp for FR-HIT
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

#include "frhit.h"

int main(int argc, char *argv[]) {
    try { 
        if (argc == 1) usage();
        Initial_Time();
        std::cerr << "Start at:  " << Curr_Time() << std::endl;
        int noptions = mGetOptions(argc, argv);
#if defined (_OPENMP)
        if (param.ncpu) omp_set_num_threads(param.ncpu);
#endif
        finDB.open(refseqFile.c_str());
        if (!finDB) {
                std::cerr << "fatal error: failed to open ref file\n";
                exit(1);
        }
        ref.Run_ConvertBinseq(finDB);
        std::cerr << "Load in " 
                  << ref.total_num 
                  << " reference seqs, total size " 
                  << ref.sum_length 
                  << " bp. " 
                  << Cal_AllTime() 
                  << " secs passed" 
                  << std::endl;
        ref.CreateIndex();
        std::cerr << "Create refseq k-mer index table. " 
                  << Cal_AllTime() 
                  << " secs passed" 
                  << std::endl;
        RunProcess();
    } catch (const std::exception & e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

