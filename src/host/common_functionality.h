/*
Copyright (c) 2019 Marius Meyer

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef COMMON_FUNCTIONALITY_H
#define COMMON_FUNCTIONALITY_H

/* C++ standard library headers */
#include <memory>

/* Project's headers */
#include "src/host/execution.h"

/*
Short description of the program
*/
#define PROGRAM_DESCRIPTION "Implementation of the LINPACK benchmark"\
                            " proposed in the HPCC benchmark suite for FPGA"

/*
Total length of the data array used for random accesses.
The array should allocate half of the available global memory space.
Keep in mind that this also depends on DATA_TYPE.
*/
#ifndef DATA_LENGTH
#define DATA_LENGTH 67108864
#endif

/*
Number of times the execution of the benchmark will be repeated.
*/
#ifndef NTIMES
#define NTIMES 1
#endif

/*
The data type used for the random accesses.
Note that it should be big enough to address the whole data array. Moreover it
has to be the same type as in the used kernels.
The signed and unsigned form of the data type have to be given separately.
*/
#ifndef DATA_TYPE
    #define DATA_TYPE cl_float
#endif

/*
Prefix of the function name of the used kernel.
It will be used to construct the full function name for the case of replications.
The full name will be
*/
#define GEFA_KERNEL "gefa"

#define ENTRY_SPACE 13

struct ProgramSettings {
    uint numRepetitions;
    uint blockSize;
    size_t matrixSize;
    bool useMemInterleaving;
    std::string kernelFileName;
};


/**
Parses and returns program options using the cxxopts library.
Supports the following parameters:
    - file name of the FPGA kernel file (-f,--file)
    - number of repetitions (-n)
    - number of kernel replications (-r)
    - data size (-d)
    - use memory interleaving
@see https://github.com/jarro2783/cxxopts

@return program settings that are created from the given program arguments
*/
std::shared_ptr<ProgramSettings>
parseProgramParameters(int argc, char * argv[]);

/**
Gaussian elemination reference implementation without pivoting.
Can be used in exchange with kernel functions for functionality testing

@param a the matrix with size of n*n
@param n size of matrix A
@param lda row with of the matrix. must be >=n

*/
void gefa_ref_nopivot(DATA_TYPE* a, ulong n, ulong lda);

/**
Solve linear equations using its LU decomposition.
Therefore solves A*x = b by solving L*y = b and then U*x = y with A = LU
where A is a matrix of size n*n

@param a the matrix a in LU representation calculated by gefa call
@param b vector b of the given equation
@param n size of matrix A
@param lda row with of the matrix. must be >=n

*/
void gesl_ref_nopivot(DATA_TYPE* a, DATA_TYPE* b, ulong n, uint lda);

/**
Print the benchmark results to stdout

@param results the struct containing the results of the benchmark execution
@param matrixSize size of the calculated matrix
*/
void printResults(std::shared_ptr<bm_execution::ExecutionResults> results,
                  size_t matrixSize);

/**
Generate a matrix using pseudo random numbers with fixed seed.
Use the matrix to generate a vector b such that
A*x = b and x = (1,1, ...,1)

@param a pointer to the matrix
@param lda width of a row in the matrix
@param n number of rows in the matrix
@param b the generated vector that holds the described condition
@param norma the maximum value in the matrix A that can be used to calculate the residual error
*/
void matgen(DATA_TYPE* a, cl_int lda, cl_int n, DATA_TYPE* b, DATA_TYPE* norma);

/**
Multiply matrix with a vector and add it to another vector.

// TODO add docs
*/
void dmxpy (int n1, DATA_TYPE* y, int n2, int ldm, DATA_TYPE* x, DATA_TYPE* m);

double checkLINPACKresults (DATA_TYPE* b_res, cl_int lda, cl_int n);

DATA_TYPE epslon (DATA_TYPE x);

int main(int argc, char * argv[]);

#endif // SRC_HOST_COMMON_FUNCTIONALITY_H_
