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

#define DATA_TYPE float

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 8
#endif

#define GEMM_BLOCK 8

/**
Load a block from global memory

@param a_block local memory buffer to store the block in
@param a the global memory buffer of the Matrix
@param x_block x position of the block
@param y_block y position of the block
@param lda_block LDA of the matrix in number of blocks
*/
void
load_block(DATA_TYPE a_block[BLOCK_SIZE][BLOCK_SIZE],
			global DATA_TYPE* restrict a,
			uint x_block, uint y_block, uint lda_block) {

	for (int i = 0; i < BLOCK_SIZE; i++) {
		#pragma unroll GLOBAL_MEM_UNROLL
		for (int j = 0; j < BLOCK_SIZE; j++) {
			a_block[i][j] = a[(y_block * lda_block * BLOCK_SIZE + x_block) * BLOCK_SIZE + j
								+ i * lda_block * BLOCK_SIZE];
		}
	}
}

/**
Store a block to global memory

@param a_block local memory buffer to load the block from
@param a the global memory buffer of the Matrix
@param x_block x position of the block
@param y_block y position of the block
@param lda_block LDA of the matrix in number of blocks
*/
void
store_block(DATA_TYPE a_block[BLOCK_SIZE][BLOCK_SIZE],
			global DATA_TYPE* restrict a,
			uint x_block, uint y_block, uint lda_block) {

	for (int i = 0; i < BLOCK_SIZE; i++) {
		#pragma unroll GLOBAL_MEM_UNROLL
		for (int j = 0; j < BLOCK_SIZE; j++) {
			a[(y_block * lda_block * BLOCK_SIZE + x_block) * BLOCK_SIZE + j
						+ i * lda_block * BLOCK_SIZE] = a_block[i][j];
		}
	}
}

/**
Calculate

c = c +  a.dot(b)

where a,b,c are matrices of size GEMM_BLOCK.
 */
void local_gemm_8x8(const DATA_TYPE a[GEMM_BLOCK][GEMM_BLOCK],
                    const DATA_TYPE b[GEMM_BLOCK][GEMM_BLOCK],
                    DATA_TYPE c_out[GEMM_BLOCK][GEMM_BLOCK]) {

    DATA_TYPE a_block[GEMM_BLOCK][GEMM_BLOCK + 1];
    DATA_TYPE b_block[GEMM_BLOCK + 1][GEMM_BLOCK];
    DATA_TYPE c_block[GEMM_BLOCK][GEMM_BLOCK];

    // Load block of matrix A and B and init C
    #pragma unroll
    for (int y=0; y<GEMM_BLOCK; y++) {
        #pragma unroll
        for (int x=0; x<GEMM_BLOCK; x++) {
            int k = (x + y) % GEMM_BLOCK;
            a_block[y][x] = a[y][k];
            b_block[y][x] = b[k][x];
            c_block[y][x] = 0;
        }
    }

    //Calculate intermediate result for block
    #pragma unroll
    for (int i=0;i<GEMM_BLOCK; i++) {
        #pragma unroll
        for (int x=0; x<GEMM_BLOCK;x++) {
            a_block[x][GEMM_BLOCK] = a_block[x][0];
            b_block[GEMM_BLOCK][x] = b_block[0][x];
        }
        #pragma unroll
        for(int y=0; y < GEMM_BLOCK; y++) {
            #pragma unroll
            for (int x=0; x<GEMM_BLOCK;x++) {
                c_block[y][x] += a_block[y][x] * b_block[y][x];
                a_block[y][x] = a_block[y][x + 1];
                b_block[y][x] = b_block[y + 1][x];
            }
        }
    }

	#pragma unroll
    for (int y=0; y<GEMM_BLOCK; y++) {
        #pragma unroll
        for (int x=0; x<GEMM_BLOCK; x++) {
            c_out[y][x] = c_block[y][x];
        }
    }
}


/**
Searches for the index of the absoulte maximum in the column and returns it.

@param column The array containing the current column
@param current_k The current column
@returns index of the absolute maximum of the values between k and BLOCK_SIZE
*/
#define BLOCK_SIZE_LOG 5
int
argmax(const DATA_TYPE column[BLOCK_SIZE], const int current_k) {
	DATA_TYPE prepared_col[BLOCK_SIZE_LOG + 1][BLOCK_SIZE];
	DATA_TYPE prepared_col_index[BLOCK_SIZE_LOG + 1][BLOCK_SIZE];
	#pragma unroll
	for (int i=0; i < BLOCK_SIZE; i++) {
		if (i < current_k) {
			prepared_col[0][i] = 0;
		} else {
			prepared_col[0][i] = fabs(column[i]);
		}
		prepared_col_index[0][i] = i;
	}
	int remaining_vals = BLOCK_SIZE;
	#pragma unroll
	for (int stage=1; stage <= BLOCK_SIZE_LOG; stage++) {
		remaining_vals = remaining_vals >> 1;
		#pragma unroll
		for (int i=0; i < remaining_vals; i++) {
			if (prepared_col[stage - 1][i] > prepared_col[stage - 1][i + remaining_vals]) {
				prepared_col[stage][i] = prepared_col[stage - 1][i];
				prepared_col_index[stage][i] = prepared_col_index[stage - 1][i];
			} else {
				prepared_col[stage][i] = prepared_col[stage - 1][i + remaining_vals];
				prepared_col_index[stage][i] = prepared_col_index[stage - 1][i + remaining_vals];
			}
		}
	}

	return prepared_col_index[BLOCK_SIZE_LOG][0];
}



/**
Standard LU factorization on a block with fixed size

Case 1 of Zhangs description

@param a_block_in Input block that has to be LU factorized
@param a_block_out Output block to write the result
@param scale_factors scaling factors that where used to scale the columns. They can be reused in C2
@param ipvt Pivoting information for C3 and solving of the system
*/
void
lu_factorization_c1(const DATA_TYPE a_block_in[BLOCK_SIZE][BLOCK_SIZE],
					DATA_TYPE a_block_out[BLOCK_SIZE][BLOCK_SIZE],
					DATA_TYPE scale_factors[BLOCK_SIZE],
					int ipvt[BLOCK_SIZE]) {

	DATA_TYPE tmp_block_write[BLOCK_SIZE][BLOCK_SIZE];
	DATA_TYPE tmp_block_read[BLOCK_SIZE][BLOCK_SIZE];

	// copy columnwise
	for (int i = 0; i < BLOCK_SIZE; i++) {
		#pragma unroll
		for (int j = 0; j <  BLOCK_SIZE; j++) {
			tmp_block_read[i][j] = a_block_in[i][j];
		}
	}

	#pragma unroll
	for (int i = 0; i < BLOCK_SIZE; i++) {
		ipvt[i] = i;
	}

	// For each diagnonal element
	#pragma max_concurrency 1
	for (int k = 0; k < BLOCK_SIZE; k++) {

		DATA_TYPE tmp_scale_col[BLOCK_SIZE];
		int col_order[BLOCK_SIZE];
		#pragma unroll
		for (int i=0; i < BLOCK_SIZE; i++) {
			col_order[i] = i;
		}

		DATA_TYPE current_col[BLOCK_SIZE];
		#pragma unroll
		for (int i=0; i<BLOCK_SIZE; i++) {
			current_col[i] = tmp_block_read[i][k];
		}
		int pivot_col = argmax(current_col, k);
		ipvt[k] = pivot_col;
		col_order[pivot_col] = k;
		col_order[k] = pivot_col;


		scale_factors[k] = -1.0 / tmp_block_read[col_order[k]][k];
		#pragma unroll
		for (int i = k + 1; i < BLOCK_SIZE; i++) {
			tmp_scale_col[i] =  current_col[col_order[i]] * scale_factors[k];
			tmp_block_write[i][k] = tmp_scale_col[i];
		}
		#pragma unroll
		for (int i = k; i < BLOCK_SIZE; i++) {
			tmp_block_write[k][i] = tmp_block_read[col_order[k]][i];
		}

		// For each column right of current diagonal element
		for (int j = k + 1; j < BLOCK_SIZE; j++) {
			// For each element below it
			#pragma unroll BLOCK_SIZE
			for (int i = 0; i < BLOCK_SIZE; i++) {
				if (i > k) {
					tmp_block_write[j][i] = tmp_block_read[col_order[j]][i] + tmp_scale_col[j] * tmp_block_read[col_order[k]][i];
				}
			}
		}
		#pragma unroll
		for (int i = k; i < BLOCK_SIZE; i++) {
			#pragma unroll
			for (int j = 0; j <  BLOCK_SIZE; j++) {
				tmp_block_read[i][j] = tmp_block_write[i][j];
			}
		}
#ifdef DEBUG
		printf("A(%d):\n", k);
		for (int j = 0; j < BLOCK_SIZE; j++) {
			for (int i = 0; i <  BLOCK_SIZE; i++) {
				printf("%f, ", tmp_block_read[j][i]);
			}
			printf("\n");
		}
		printf("\n\n");
#endif
	}
	for (int i = 0; i < BLOCK_SIZE; i++) {
		#pragma unroll
		for (int j = 0; j <  BLOCK_SIZE; j++) {
			a_block_out[i][j] = tmp_block_read[i][j];
		}
	}

}

/**
Modifying the blocks on the leftmost side

Case 2 of Zhangs description

@param top_block LU factorized top block
@param current_block_in Current input block
@param current_block_out Block to write the output to
@param scale_factors Scale factors that were calculated during LU factorization
*/
void
left_blocks_c2(const DATA_TYPE top_block[BLOCK_SIZE][BLOCK_SIZE],
				const DATA_TYPE current_block_in[BLOCK_SIZE][BLOCK_SIZE],
				DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE],
				const DATA_TYPE scale_factors[BLOCK_SIZE]) {

	DATA_TYPE tmp_block_write2[BLOCK_SIZE][BLOCK_SIZE];
	DATA_TYPE tmp_block_read2[BLOCK_SIZE][BLOCK_SIZE];
	DATA_TYPE tmp_scale_col[BLOCK_SIZE];

	// copy columnwise
	for (int i = 0; i < BLOCK_SIZE; i++) {
		#pragma unroll
		for (int j = 0; j <  BLOCK_SIZE; j++) {
			tmp_block_read2[i][j] = current_block_in[j][i];
		}
	}
	// For each diagonal element in top block
	#pragma max_concurrency 1
	for (int k=0; k < BLOCK_SIZE; k++) {
		// For each element below it in current block
		#pragma unroll
		for (int i=0; i < BLOCK_SIZE; i++) {
			// printf("C2: %f * %f\n",tmp_block2[i][k], scale_factors[k]);
			tmp_scale_col[i] = tmp_block_read2[k][i] * scale_factors[k];
			tmp_block_write2[k][i] = tmp_scale_col[i];
		}
		// For each column right of the current diagnonal element
		for (int j = k+1; j < BLOCK_SIZE; j++) {
			#pragma unroll
			for (int i = 0; i < BLOCK_SIZE; i++) {
				tmp_block_write2[j][i] =
							tmp_block_read2[j][i] + tmp_scale_col[i]
												* top_block[k][j];
			}
		}
		for (int i = k; i < BLOCK_SIZE; i++) {
			#pragma unroll
			for (int j = 0; j <  BLOCK_SIZE; j++) {
				tmp_block_read2[i][j] = tmp_block_write2[i][j];
			}
		}
	}
	for (int i = 0; i < BLOCK_SIZE; i++) {
		#pragma unroll
		for (int j = 0; j <  BLOCK_SIZE; j++) {
			current_block_out[j][i] = tmp_block_write2[i][j];
		}
	}
}

/**
Modifying the blocks on the top but not on the left

Case 3 of Zhangs description

@param left_block LU factorized left block
@param current_block_in Current input block
@param current_block_out Block to write the output to
@param ipvt Pivot information created by the LU factorization
*/
void
top_blocks_c3(const DATA_TYPE left_block[BLOCK_SIZE][BLOCK_SIZE],
			  const DATA_TYPE current_block_in[BLOCK_SIZE][BLOCK_SIZE],
			  DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE],
			  int ipvt[BLOCK_SIZE]) {
	DATA_TYPE tmp_block_read3[BLOCK_SIZE][BLOCK_SIZE];
	DATA_TYPE tmp_block_write3[BLOCK_SIZE][BLOCK_SIZE];


	for (int j = 0; j < BLOCK_SIZE; j++) {
		#pragma unroll
		for (int i = 0; i <  BLOCK_SIZE; i++) {
			tmp_block_read3[j][i] = current_block_in[j][i];
		}
	}

	// For each diagonal element in left block
	#pragma max_concurrency 1
	for (int k=0; k < BLOCK_SIZE; k++) {
		uint col_order[BLOCK_SIZE];
		#pragma unroll
		for (int i=0; i < BLOCK_SIZE; i++) {
			col_order[i] = i;
		}
		col_order[k] = ipvt[k];
		col_order[ipvt[k]] = k;
		// For each column in current block
		for (int j = k; j < BLOCK_SIZE; j++) {
			DATA_TYPE multiply = 0.0;
			// For each element below it
			if (j > k) {
				multiply = left_block[j][k];
			}
			#pragma unroll
			for (int i = 0; i < BLOCK_SIZE; i++) {
				tmp_block_write3[j][i] = tmp_block_read3[col_order[j]][i]
					+ multiply * tmp_block_read3[col_order[k]][i];
			}
		}
		for (int j = k; j < BLOCK_SIZE; j++) {
			#pragma unroll
			for (int i = 0; i <  BLOCK_SIZE; i++) {
				tmp_block_read3[j][i] = tmp_block_write3[j][i];
			}
		}
	}
	for (int j = 0; j < BLOCK_SIZE; j++) {
		#pragma unroll
		for (int i = 0; i <  BLOCK_SIZE; i++) {
			current_block_out[j][i] = tmp_block_write3[j][i];
		}
	}
}

/**
Modifying the inner blocks

Case 4 of Zhangs description

@param left_block Most left block that was modified by C2 before
@param top_block Most upper block that was modified by C3 before
@param current_block_in Current input block
@param current_block_out Block to write the output to
*/
void
inner_blocks_c4(const DATA_TYPE left_block[BLOCK_SIZE][BLOCK_SIZE],
				const DATA_TYPE top_block[BLOCK_SIZE][BLOCK_SIZE],
				const DATA_TYPE current_block_in[BLOCK_SIZE][BLOCK_SIZE],
				DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE]) {
	DATA_TYPE tmp_top_block[BLOCK_SIZE / GEMM_BLOCK][BLOCK_SIZE / GEMM_BLOCK]
							 [GEMM_BLOCK][GEMM_BLOCK];
	DATA_TYPE tmp_left_block[BLOCK_SIZE / GEMM_BLOCK][BLOCK_SIZE / GEMM_BLOCK]
							 [GEMM_BLOCK][GEMM_BLOCK];

	#pragma loop_coalesce
	for (int i = 0; i < BLOCK_SIZE / GEMM_BLOCK; i++) {
		for (int j = 0; j < BLOCK_SIZE / GEMM_BLOCK; j++) {
			#pragma unroll
			for (int ii = 0; ii < GEMM_BLOCK; ii++) {
				#pragma unroll
				for (int jj = 0; jj < GEMM_BLOCK; jj++) {
					tmp_top_block[i][j][ii][jj] = top_block[i * GEMM_BLOCK + ii]
														[j * GEMM_BLOCK + jj];
					tmp_left_block[i][j][ii][jj] = left_block[i * GEMM_BLOCK + ii]
														[j * GEMM_BLOCK + jj];
				}
			}
		}
	}

	#pragma loop_coalesce 3
	// For each diagonal element in left block
	for (int k=0; k < BLOCK_SIZE / GEMM_BLOCK; k++) {
		// For each column in top block
		for (int i = 0; i < BLOCK_SIZE / GEMM_BLOCK; i++) {
			// For each element below it in current block
			for (int j = 0; j < BLOCK_SIZE / GEMM_BLOCK; j++) {
				DATA_TYPE   tmp_block_out[GEMM_BLOCK][GEMM_BLOCK];
				local_gemm_8x8(tmp_left_block[i][k], tmp_top_block[k][j], tmp_block_out);
				if (k == 0) {
					for (int ii = 0; ii < GEMM_BLOCK; ii++) {
						#pragma unroll
						for (int jj = 0; jj < GEMM_BLOCK; jj++) {
							current_block_out[i * GEMM_BLOCK + ii][j * GEMM_BLOCK + jj] = current_block_in[i * GEMM_BLOCK + ii][j * GEMM_BLOCK + jj] + tmp_block_out[ii][jj];
						}
					}
				} else {
					for (int ii = 0; ii < GEMM_BLOCK; ii++) {
						#pragma unroll
						for (int jj = 0; jj < GEMM_BLOCK; jj++) {
							current_block_out[i * GEMM_BLOCK + ii][j * GEMM_BLOCK + jj] += tmp_block_out[ii][jj];
						}
					}
				}
			}
		}
	}
}


/**
LU factorization kernel

@param a The data array representing the whole matrix in global memory
@param pvt Pivoting information
@param a_size the x and y size of the matrix in blocks
*/
__attribute__((uses_global_work_offset(0)))
__kernel
void gefa(global DATA_TYPE* restrict a, global int* restrict pvt,  uint a_size) {

	// For each diagonal block do the following
	#pragma disable_loop_pipelining
	for (int diagonal_block=0; diagonal_block < a_size; diagonal_block++) {
		DATA_TYPE diag_block[BLOCK_SIZE][BLOCK_SIZE];
		DATA_TYPE diag_block_out[BLOCK_SIZE][BLOCK_SIZE];
		// load next block for factorization
		load_block(diag_block, a, diagonal_block, diagonal_block, a_size);

		DATA_TYPE scale_factors[BLOCK_SIZE];
		int ipvt[BLOCK_SIZE];

		// execute factorization of next block
		lu_factorization_c1(diag_block, diag_block_out, scale_factors,
													ipvt);

		#pragma unroll GLOBAL_MEM_UNROLL
		for (int i=0; i<BLOCK_SIZE; i++) {
			pvt[diagonal_block * BLOCK_SIZE + i] = diagonal_block * BLOCK_SIZE
																+ ipvt[i];
		}

		store_block(diag_block_out, a, diagonal_block, diagonal_block, a_size);

		for (int inner_block = diagonal_block + 1; inner_block < a_size;
			inner_block++) {
			DATA_TYPE top_block[BLOCK_SIZE][BLOCK_SIZE];
			DATA_TYPE left_block[BLOCK_SIZE][BLOCK_SIZE];
			DATA_TYPE top_block_out[BLOCK_SIZE][BLOCK_SIZE];
			DATA_TYPE left_block_out[BLOCK_SIZE][BLOCK_SIZE];
			// update top block
			load_block(left_block, a, diagonal_block,
										inner_block, a_size);
			load_block(top_block, a, inner_block, diagonal_block, a_size);
			left_blocks_c2(diag_block_out, left_block,
								left_block_out, scale_factors);
			top_blocks_c3(diag_block_out, top_block, top_block_out, ipvt);
			store_block(top_block_out, a, inner_block,
													diagonal_block, a_size);
			store_block(left_block_out, a, diagonal_block,
										inner_block, a_size);

		}

		// process all blocks columnwise
		for (int inner_x_block = diagonal_block + 1; inner_x_block < a_size;
			inner_x_block++) {
			DATA_TYPE top_block_out[BLOCK_SIZE][BLOCK_SIZE];
			// update top block
			load_block(top_block_out, a, inner_x_block, diagonal_block, a_size);

			for (int inner_y_block = diagonal_block + 1;
								inner_y_block < a_size; inner_y_block++) {
				DATA_TYPE left_block_out[BLOCK_SIZE][BLOCK_SIZE];
				DATA_TYPE current_block[BLOCK_SIZE][BLOCK_SIZE];
				DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE];

				load_block(left_block_out, a, diagonal_block,
											inner_y_block, a_size);

				// update inner block
				load_block(current_block, a, inner_x_block,
												inner_y_block, a_size);

				inner_blocks_c4(left_block_out, top_block_out, current_block,
													current_block_out);

				store_block(current_block_out, a, inner_x_block,
												inner_y_block, a_size);
			}
		}
	}
}
