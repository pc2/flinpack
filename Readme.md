# LINPACK for FPGA

This repository contains the LINPACK for FPGA and its OpenCL kernels.
Currently only the  Intel FPGA SDK for OpenCL utility is supported.

The implementation is currently work in progess and is not feature complete.
Read the section **Implementation Details** for more information.


## Build

The code has the following dependencies:

- C++ compiler (GCC >= 4.9 with libc++)
- Intel FPGA SDK for OpenCL (tested with 19.2)

Depending on the use case you may want to change certain lines in the
Makefile:

1. Check the location of the used compilers (A C++ compiler and the aoc/aocl)
2. Update the board name in the Makefile or give the new board name as an
   argument `BOARD` to make.
3. Set the block size with `BLOCK_SIZE`. If modified you should also give the
   logarithm of the size in `BLOCK_SIZE_LOG`. It is used in the kernel code
   for the maximum determination.
4. Build the host program or the kernels by using the available build targets.

For more detailed information about the available build targets call:

    make

without specifying a target.
This will print a list of the targets together with a short description.

To compile all necessary bitstreams and executables for the benchmark run:

    make host kernel

To make it easier to generate different versions of the kernels, it
is possible to specify a variable `BUILD_SUFFIX` when executing make.
This suffix will be added to the kernel name after generation.

Example:

	make host BUILD_SUFFIX=19.2

Will build the host and name the binary after the given build suffix.
So the host would be named `execution_blocked_pvt_19.2`.
The file will be placed in a folder called `bin` in the root of this project.

Moreover it is possible to specifiy additional aoc compiler parameters using the
`AOC_FLAGS` variable.
For example:

    make kernel AOC_FLAGS="-fpc -fp-relaxed"

## Execution

The created host needs the kernel file as a program argument.
The kernel file has to be specified using the `-f` flag like the following:

    ./execution_blocked_pvt_19.2 -f path/to/file.aocx

It is also possible to give additional settings. To get a more detailed overview
of the available settings execute:

    ./execution_blocked_pvt_19.2 -h

## Implementation Details

The benchmark will measure the elapsed time to execute a kernel for performing
an LU factorization.
It will use the time to calculate the FLOP/s.
Buffer transfer is currently not measured.
The solving of the linear equations is currently done on the CPU.

The updates are done unaligned and randomly directly on the global memory.
The repository contains two different implementations:
- `blocked`: A blocked, unoptimized kernel that performs the LU factorization
   without pivoting.
- `blocked_pvt`: Blocked kernel that performs the LU factorization with pivoting
   over the whole block.

#### Adjustable Parameters

The following table shows the modifiable parameters and in which kernel they
can be used to adjust the code.
The parameter `MATRIX_SIZE` can also be adjusted at
runtime in the host using command line parameters.
Use the `-h` option to show available parameters when running the host code.
You can specify the kernel type using the `TYPE` parameter when running make.
See the example below the table.

| Parameter         | `blocked`/<br>`blocked_pvt`/<br>Host      | Details                                  |
|------------------ | ------------------------------------------------------ | ---------------------------------------- |
| `TYPE`           |:white_check_mark:/:white_check_mark:/:white_check_mark:     |   Type of the used kernel. Default is `blocked_pvt`.  |
| `BOARD`           |:white_check_mark:/:white_check_mark:/:x:     |   Name of the target board               |
| `BUILD_SUFFIX`    |:white_check_mark:/:white_check_mark:/:white_check_mark:| Addition to the kernel name              |
| `AOC_FLAGS`       |:white_check_mark:/:white_check_mark:/:x:               | Additional compile flags for `aoc`       |
| `MATRIX_SIZE` |:x:/:x:/:white_check_mark:                              | Default matrix size. Can also be specified in the host at runtime.   |
| `BLOCK_SIZE`    |:white_check_mark:/:white_check_mark:/:white_check_mark:             | Size of a block.  |
| `BLOCK_SIZE_LOG`    |:x:/:white_check_mark:/:x:             | Log2 of the size of a block.  |
| `GLOBAL_MEM_UNROLL`|:white_check_mark:/:white_check_mark:/:x:              | Unrolling of loops that access the global memory |
| `CXX_FLAGS`       |:x:/:x:/:white_check_mark:                              | Additional C++ compiler flags            |

Example for synthesizing a kernel to create a profiling report:

```bash
make kernel_profile BOARD=p520_hpc_sg280l CXX_FLAGS=-03 \
AOC_FLAGS="-fpc -fp-relaxed" BLOCK_SIZE=32 BLOCK_SIZE_LOG=5 TYPE=blocked_pvt
```

#### Work in Progress

The implementation is currently work in progress and currently only covers the
GEFA calculation on FPGA.
A rough overview of the WIP with focus on the pivoting kernel:

- Routines C1 to C3 are not optimized and C4 reduces fMax.
- Only block-wise partial pivoting is used instead of partial pivoting over
  the whole matrix. This increases the error in the calculation.
- GESL not implemented on FPGA.


## Result Interpretation

The host code will print the results of the execution to the standard output.
The result  summary looks similar to this:

    norm. resid        resid       machep       x[0]-1     x[n-1]-1
    9.40193e+00  2.87056e-04  1.19209e-07 -3.21865e-05  2.57969e-04
    best         mean         GFLOPS       error
    1.57262e-01  1.57262e-01  7.11221e-02  9.40193e+00

The first row contains data from the correctness check that is done once when
executing the benchmark:
- `resid`: The maximum residual error when multiplying the result vector with
   the matrix and subtract by the expected result.
- `norm. resid`: The normalized residual error based on `resid`.
- `machep`: machine epsilon that gives an upper bound for rounding errors due
   to the used floating point format.
- `x[0] - 1`: The first element of the result vector minus 1. It should be
   close to 0. The same holds for `x[n-1] - 1` which is the last element of the
   vector.

The second row contains the measured performance of the benchmark:
- `best`: The best measured time for executing the benchmark in seconds.
- `mean`: The arithmetic mean of all measured execution times in seconds.
- `GFLOPS`: GFLOP/s achieved for the calculation using the best measured time.
- `error`: Same as `norm. resid` to complete the performance overview.
