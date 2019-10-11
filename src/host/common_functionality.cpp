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

/* Related header files */
#include "src/host/common_functionality.h"

/* C++ standard library headers */
#include <iostream>
#include <cmath>
#include <string>
#include <limits>
#include <iomanip>
#include <memory>
#include <vector>

/* External library headers */
#include "CL/cl.hpp"
#if QUARTUS_MAJOR_VERSION > 18
#include "CL/cl_ext_intelfpga.h"
#endif

/* Project's headers */
#include "src/host/benchmark_helper.h"
#include "src/host/cxxopts.hpp"
#include "src/host/execution.h"


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
parseProgramParameters(int argc, char * argv[]) {
    // Defining and parsing program options
    cxxopts::Options options(argv[0], PROGRAM_DESCRIPTION);
    options.add_options()
        ("f,file", "Kernel file name", cxxopts::value<std::string>())
        ("n", "Number of repetitions",
                cxxopts::value<uint>()->default_value(std::to_string(NTIMES)))
        ("b", "Used block size",
            cxxopts::value<uint>()->default_value(std::to_string(BLOCK_SIZE)))
        ("m,matrix", "Size of the matrix (NxN)",
                cxxopts::value<size_t>()
                                ->default_value(std::to_string(MATRIX_SIZE)))
        ("i,nointerleaving", "Disable memory interleaving")
        ("h,help", "Print this help");
    cxxopts::ParseResult result = options.parse(argc, argv);

    // Check parsed options and handle special cases
    if (result.count("f") <= 0) {
        // Path to the kernel file is mandatory - exit if not given!
        std::cerr << "Kernel file must be given! Aborting" << std::endl;
        std::cout << options.help() << std::endl;
        exit(1);
    }
    if (result.count("h")) {
        // Just print help when argument is given
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Create program settings from program arguments
    std::shared_ptr<ProgramSettings> sharedSettings(
            new ProgramSettings {result["n"].as<uint>(), result["b"].as<uint>(),
                                result["m"].as<size_t>(),
                                static_cast<bool>(result.count("i") <= 0),
                                result["f"].as<std::string>()});
    return sharedSettings;
}

/**
Print the benchmark Results

@param results The result struct provided by the calculation call
@param dataSize The size of the used data array

*/
void printResults(std::shared_ptr<bm_execution::ExecutionResults> results,
                  size_t dataSize) {
    std::cout << std::setw(ENTRY_SPACE)
              << "best" << std::setw(ENTRY_SPACE) << "mean"
              << std::setw(ENTRY_SPACE) << "GUOPS"
              << std::setw(ENTRY_SPACE) << "error" << std::endl;

    // Calculate performance for kernel execution plus data transfer
    double tmean = 0;
    double tmin = std::numeric_limits<double>::max();
    double gops = ((2.0e0*(dataSize*dataSize*dataSize))/3.0
                    + 2.0*(dataSize*dataSize)) / 1.0e9;
    for (double currentTime : results->times) {
        tmean +=  currentTime;
        if (currentTime < tmin) {
            tmin = currentTime;
        }
    }
    tmean = tmean / results->times.size();

    std::cout << std::setw(ENTRY_SPACE)
              << tmin << std::setw(ENTRY_SPACE) << tmean
              << std::setw(ENTRY_SPACE) << gops / tmin
              << std::setw(ENTRY_SPACE) << (results->errorRate)
              << std::endl;
}

void matgen(DATA_TYPE* a, cl_int lda, cl_int n, DATA_TYPE* b,
            DATA_TYPE* norma) {
    int init = 1325;
    *norma = 0.0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            init = 3125*init % 65536;
            a[lda*j+i] = (init - 32768.0)/16384.0;
            *norma = (a[lda*j+i] > *norma) ? a[lda*j+i] : *norma;
        }
        for (int i = n; i < lda; i++) {
            a[lda*j+i] = 0;
        }
    }
    for (int i = 0; i < n; i++) {
          b[i] = 0.0;
      }
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            b[j] += a[lda*j+i];
        }
    }
}

/**
Standard LU factorization on a block with fixed size

Case 1 of Zhangs description
*/
void
gefa_ref_nopivot(DATA_TYPE* a, ulong n, ulong lda) {
    // For each diagnonal element
    for (int k = 0; k < n - 1; k++) {
        // For each element below it
        for (int i = k + 1; i < n; i++) {
            a[i * lda + k] *= 1.0 / a[k * lda + k];
        }
        // For each column right of current diagonal element
        for (int j = k + 1; j < n; j++) {
            // For each element below it
            for (int i = k+1; i < n; i++) {
                a[i * lda + j] -= a[i * lda + k] * a[k * lda + j];
            }
        }
    }
}

void
gesl_ref_nopivot(DATA_TYPE* a, DATA_TYPE* b, ulong n, uint lda) {
    // solve l*y = b
    // For each row in matrix
    for (int k = 0; k < n-1; k++) {
        // For each row below add
        for (int i = k+1; i < n; i++) {
            // add solved upper row to current row
            b[i] -= b[k] * a[lda*i + k];
        }
    }

    // now solve  u*x = y

    for (int k = n-1; k >= 0; k--) {
        b[k] = b[k]/a[lda*k + k];
        DATA_TYPE t = -b[k];
        for (int i = 0; i < k; i++) {
            b[i] += t * a[lda*i + k];
        }
    }
}

void dmxpy(int n1, DATA_TYPE* y, int n2, int ldm, DATA_TYPE* x, DATA_TYPE* m) {
    int j, i, jmin;
    /* cleanup odd vector */

    j = n2 % 2;
    if (j >= 1) {
        j = j - 1;
        for (i = 0; i < n1; i++)
            y[i] = (y[i]) + x[j]*m[ldm*j+i];
    }

    /* cleanup odd group of two vectors */

    j = n2 % 4;
    if (j >= 2) {
        j = j - 1;
        for (i = 0; i < n1; i++)
            y[i] = ( (y[i]) + x[j-1]*m[ldm*(j-1)+i]) + x[j]*m[ldm*j+i];
    }

    /* cleanup odd group of four vectors */

    j = n2 % 8;
    if (j >= 4) {
        j = j - 1;
        for (i = 0; i < n1; i++)
            y[i] = ((( (y[i])
                    + x[j-3]*m[ldm*(j-3)+i])
                    + x[j-2]*m[ldm*(j-2)+i])
                    + x[j-1]*m[ldm*(j-1)+i]) + x[j]*m[ldm*j+i];
    }

    /* cleanup odd group of eight vectors */

    j = n2 % 16;
    if (j >= 8) {
        j = j - 1;
        for (i = 0; i < n1; i++)
        y[i] = ((((((( (y[i])
                + x[j-7]*m[ldm*(j-7)+i]) + x[j-6]*m[ldm*(j-6)+i])
                + x[j-5]*m[ldm*(j-5)+i]) + x[j-4]*m[ldm*(j-4)+i])
                + x[j-3]*m[ldm*(j-3)+i]) + x[j-2]*m[ldm*(j-2)+i])
                + x[j-1]*m[ldm*(j-1)+i]) + x[j]  *m[ldm*j+i];
    }

    /* main loop - groups of sixteen vectors */

    jmin = (n2%16)+16;
    for (j = jmin-1; j < n2; j = j + 16) {
        for (i = 0; i < n1; i++)
            y[i] = ((((((((((((((( (y[i])
                    + x[j-15]*m[ldm*(j-15)+i])
                    + x[j-14]*m[ldm*(j-14)+i])
                    + x[j-13]*m[ldm*(j-13)+i])
                    + x[j-12]*m[ldm*(j-12)+i])
                    + x[j-11]*m[ldm*(j-11)+i])
                    + x[j-10]*m[ldm*(j-10)+i])
                    + x[j- 9]*m[ldm*(j- 9)+i])
                    + x[j- 8]*m[ldm*(j- 8)+i])
                    + x[j- 7]*m[ldm*(j- 7)+i])
                    + x[j- 6]*m[ldm*(j- 6)+i])
                    + x[j- 5]*m[ldm*(j- 5)+i])
                    + x[j- 4]*m[ldm*(j- 4)+i])
                    + x[j- 3]*m[ldm*(j- 3)+i])
                    + x[j- 2]*m[ldm*(j- 2)+i])
                    + x[j- 1]*m[ldm*(j- 1)+i])
                    + x[j]   *m[ldm*j+i];
    }
}

double
checkLINPACKresults(DATA_TYPE* b_res, cl_int lda, cl_int n) {
    DATA_TYPE a[lda*n];
    DATA_TYPE norma = 0;
    DATA_TYPE x[n];
    DATA_TYPE b[n];
    /*     compute a residual to verify results.  */

    for (int i = 0; i < n; i++) {
        x[i] = b_res[i];
        b[i] = b_res[i];
    }
    matgen(a, lda, n, b, &norma);
    for (int i = 0; i < n; i++) {
        b[i] = -b[i];
    }
    dmxpy(n, b, n, lda, x, a);
    DATA_TYPE resid = 0.0;
    DATA_TYPE normx = 0.0;
    for (int i = 0; i < n; i++) {
        resid = (resid > fabs(b[i])) ? resid : fabs((DATA_TYPE)b[i]);
        normx = (normx > fabs(x[i])) ? normx : fabs(x[i]);
    }

    DATA_TYPE eps = epslon(1.0);
    DATA_TYPE residn = resid / (n*norma*normx*eps);

    std::cout << "  norm. resid        resid       "\
                 "machep       x[0]-1     x[n-1]-1" << std::endl;
    std::cout << std::setw(ENTRY_SPACE) << residn << std::setw(ENTRY_SPACE)
              << resid << std::setw(ENTRY_SPACE) << eps
              << std::setw(ENTRY_SPACE) << x[0]-1 << std::setw(ENTRY_SPACE)
              << x[n-1]-1 << std::endl;
    return residn;
}

DATA_TYPE epslon(DATA_TYPE x) {
    DATA_TYPE a, b, c, eps;

    a = 4.0e0/3.0e0;
    eps = 0.0;
    while (eps == 0.0) {
        b = a - 1.0;
        c = b + b + b;
        eps = fabs(static_cast<double>(c-1.0));
    }
    return (eps*fabs(static_cast<double>(x)));
}


/**
The program entry point.
Prepares the FPGA and executes the kernels on the device.
*/
int main(int argc, char * argv[]) {
    // Setup benchmark
    std::shared_ptr<ProgramSettings> programSettings =
                                            parseProgramParameters(argc, argv);
    bm_helper::setupEnvironmentAndClocks();
    std::vector<cl::Device> usedDevice = bm_helper::selectFPGADevice();
    cl::Context context = cl::Context(usedDevice);
    const char* usedKernel = programSettings->kernelFileName.c_str();
    cl::Program program = bm_helper::fpgaSetup(context, usedDevice, usedKernel);

    // Give setup summary
    std::cout << "Summary:" << std::endl
              << "Kernel Repetitions:  " << programSettings->numRepetitions
              << std::endl
              << "Block size:          " << programSettings->blockSize
              << std::endl
              << "Total matrix size:   " << programSettings->matrixSize
              << std::endl
              << "Memory Interleaving: " << programSettings->useMemInterleaving
              << std::endl
              << "Kernel file:         " << programSettings->kernelFileName
              << std::endl
              << "Device:              "
              << usedDevice[0].getInfo<CL_DEVICE_NAME>() << std::endl
              << HLINE
              << "Start benchmark using the given configuration." << std::endl
              << HLINE;

    // Start actual benchmark
    auto results = bm_execution::calculate(context, usedDevice[0], program,
              programSettings->numRepetitions, programSettings->matrixSize);

    printResults(results, programSettings->matrixSize);

    return 0;
}
