##
#  Author: Marius Meyer
#  eMail:  marius.meyer@upb.de
#  Date:   2019/07/24
#
#  This makefile compiles the Random Access benchmark for FPGA and its OpenCL kernels.
#  Currently only the  Intel(R) FPGA SDK for OpenCL(TM) utitlity is supported.
#
#  Please read the README contained in this folder for more information and
#  instructions how to compile and execute the benchmark.

# Used compilers for C code and OpenCL kernels
CXX := g++
AOC := aoc
AOCL := aocl
MKDIR_P := mkdir -p


$(info ***************************)
$(info Collected information)

# Check Quartus version
ifndef QUARTUS_VERSION
$(info QUARTUS_VERSION not defined! Quartus is not set up correctly or version is too old (<17.1). In the latter case:)
$(info Define the variable with the used Quartus version by giving it as an argument to make: e.g. make QUARTUS_VERSION=17.1.0)
$(error QUARTUS_VERSION not defined!)
else
$(info QUARTUS_VERSION         = $(QUARTUS_VERSION))
endif

QUARTUS_MAJOR_VERSION := $(shell echo $(QUARTUS_VERSION) | cut -d "." -f 1)

# OpenCL compile and link flags.
AOCL_COMPILE_CONFIG := $(shell $(AOCL) compile-config )
AOCL_LINK_CONFIG := $(shell $(AOCL) link-config )

BIN_DIR := bin/
SRC_DIR := src/

ifdef BUILD_SUFFIX
	EXT_BUILD_SUFFIX := _$(BUILD_SUFFIX)
endif

## Build settings
#
# This section contains all build parameters that might be considered to
# change for compilation and synthesis.
#
##
TYPE := blocked
GLOBAL_MEM_UNROLL := 16
BOARD := p520_hpc_sg280l
MATRIX_SIZE := 256
BLOCK_SIZE := 16
## End build settings

# The source files that differ between the chosen type
MAIN_SRC := execution_$(TYPE).cpp
KERNEL_MAIN_SRC := lu_$(TYPE).cl

KERNEL_SRC := $(SRC_DIR)device/$(KERNEL_MAIN_SRC)
SRCS := $(patsubst %, $(SRC_DIR)host/%, $(MAIN_SRC) benchmark_helper.cpp common_functionality.cpp)
TARGET := $(MAIN_SRC:.cpp=)$(EXT_BUILD_SUFFIX)
KERNEL_TARGET := $(KERNEL_MAIN_SRC:.cl=)$(EXT_BUILD_SUFFIX)

COMMON_FLAGS := -DBLOCK_SIZE=$(BLOCK_SIZE)\
 				-DQUARTUS_MAJOR_VERSION=$(QUARTUS_MAJOR_VERSION)
CXX_PARAMS := $(CXX_FLAGS) -DMATRIX_SIZE=$(MATRIX_SIZE)
AOC_PARAMS := $(AOC_FLAGS) -board=$(BOARD) -DGLOBAL_MEM_UNROLL=$(GLOBAL_MEM_UNROLL)

ifdef DATA_TYPE
CXX_PARAMS += -DDATA_TYPE=cl_$(DATA_TYPE) -DDATA_TYPE_UNSIGNED=cl_$(DATA_TYPE_UNSIGNED)
AOC_PARAMS += -DDATA_TYPE=$(DATA_TYPE) -DDATA_TYPE_UNSIGNED=$(DATA_TYPE_UNSIGNED)
endif

CXX_PARAMS += -I. -I./cxxopts/include --std=c++11

$(info Common Parameters:)
$(info BUILD_SUFFIX            = $(BUILD_SUFFIX))
$(info MATRIX_SIZE             = $(MATRIX_SIZE))
$(info BLOCK_SIZE              = $(BLOCK_SIZE))
$(info TYPE                    = $(TYPE))
$(info Device Only Parameters:)
$(info BOARD                   = $(BOARD))
$(info AOC_FLAGS               = $(AOC_FLAGS))
$(info GLOBAL_MEM_UNROLL       = $(GLOBAL_MEM_UNROLL))
$(info Host Only Parameters:)
$(info CXX_FLAGS               = $(CXX_FLAGS))
$(info ***************************)

default: info
	$(error No target specified)

info:
	$(info *************************************************)
	$(info Please specify one ore more of the listed targets)
	$(info *************************************************)
	$(info Host Code:)
	$(info host                         = Use memory interleaving to store the arrays on the FPGA)
	$(info *************************************************)
	$(info Kernels:)
	$(info kernel                       = Compile global memory kernel)
	$(info kernel_emulate               = Compile  global memory kernel for emulation)
	$(info kernel_profile               = Compile  global memory kernel with profiling information enabled)
	$(info run_emu                      = Creates host and kernel_emulate and executes the emulation with GDB)
	$(info ************************************************)
	$(info info                         = Print this list of available targets)
	$(info ************************************************)
	$(info Additional compile flags for the kernels can be provided in AOC_FLAGS.)
	$(info To disable memory interleaving: make kernel AOC_FLAGS=-no-interleaving=default)


host: $(SRCS)
	$(MKDIR_P) $(BIN_DIR)
	$(CXX) $(CXX_PARAMS) $(AOCL_COMPILE_CONFIG) $(COMMON_FLAGS)\
	$(SRCS) $(AOCL_LINK_CONFIG) -o $(BIN_DIR)$(TARGET)

kernel: $(KERNEL_SRC)
	$(MKDIR_P) $(BIN_DIR)
	$(AOC) $(AOC_PARAMS) $(COMMON_FLAGS) -o $(BIN_DIR)$(KERNEL_TARGET) $(KERNEL_SRC)

kernel_emulate: $(KERNEL_SRC)
	$(MKDIR_P) $(BIN_DIR)
	$(AOC) -march=emulator $(AOC_PARAMS) $(COMMON_FLAGS) -o $(BIN_DIR)$(KERNEL_TARGET)_emulate $(KERNEL_SRC)


run_emu: host kernel_emulate
	chmod +x $(BIN_DIR)$(TARGET)
	CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA=1 gdb --args $(BIN_DIR)$(TARGET) -f $(BIN_DIR)$(KERNEL_TARGET)_emulate.aocx

kernel_profile: $(KERNEL_SRC)
	$(MKDIR_P) $(BIN_DIR)
	$(AOC) $(AOC_PARAMS) $(COMMON_FLAGS) -profile -o $(BIN_DIR)$(KERNEL_TARGET)_profile $(KERNEL_SRC)

kernel_report: $(KERNEL_SRC)
	$(MKDIR_P) $(BIN_DIR)
ifeq ($(QUARTUS_MAJOR_VERSION), 17)
	$(AOC) $(AOC_PARAMS) $(COMMON_FLAGS) -c -o $(BIN_DIR)$(KERNEL_TARGET)_report $(KERNEL_SRC) -report
else
	$(AOC) $(AOC_PARAMS) $(COMMON_FLAGS) -rtl -o $(BIN_DIR)$(KERNEL_TARGET)_report $(KERNEL_SRC) -report
endif

cleanhost:
	rm -f $(BIN_DIR)$(TARGET)

cleanall: cleanhost
	rm -rf $(BIN_DIR)
