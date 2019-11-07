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

/**
Load a block from global memory

@param a_block local memory buffer to store the block in
@param a the global memory buffer of the Matrix
@param x_block x position of the block
@param y_block y position of the block
@param lda_block LDA of the matrix in number of blocks
*/
void
load_block(local DATA_TYPE a_block[BLOCK_SIZE][BLOCK_SIZE],
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
store_block(local DATA_TYPE a_block[BLOCK_SIZE][BLOCK_SIZE],
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
Copy a block to another memory location

@param a_from local memory buffer to load the block from
@param a_to the local memory buffer to copy the block into
*/
void
copy_block(local const DATA_TYPE a_from[BLOCK_SIZE][BLOCK_SIZE],
			DATA_TYPE a_to[BLOCK_SIZE][BLOCK_SIZE]) {

	for (int i = 0; i < BLOCK_SIZE; i++) {
		#pragma unroll
		for (int j = 0; j < BLOCK_SIZE; j++) {
			a_to[i][j] = a_from[i][j];
		}
	}
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
lu_factorization_c1(local const DATA_TYPE a_block_in[BLOCK_SIZE][BLOCK_SIZE],
					local DATA_TYPE a_block_out[BLOCK_SIZE][BLOCK_SIZE],
					DATA_TYPE scale_factors[BLOCK_SIZE],
					uint ipvt[BLOCK_SIZE]) {

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

		for (int i=0; i < BLOCK_SIZE; i++) {
			col_order[i] = i;
		}

		DATA_TYPE max_val = fabs(tmp_block_read[k][k]);
		int pivot_col = k;
		for (int i=k+1; i < BLOCK_SIZE; i++) {
			DATA_TYPE curr_val = fabs(tmp_block_read[i][k]);
			if (curr_val > max_val) {
				max_val = curr_val;
				pivot_col = i;
			}
		}
		ipvt[k] = pivot_col;
		col_order[pivot_col] = k;
		col_order[k] = pivot_col;


		scale_factors[k] = 1.0 / tmp_block_read[col_order[k]][k];
		#pragma unroll
		for (int i = k + 1; i < BLOCK_SIZE; i++) {
			tmp_scale_col[i] =  tmp_block_read[col_order[i]][k] * scale_factors[k];
			tmp_block_write[i][k] = tmp_scale_col[i];
		}
		#pragma unroll
		for (int i = k + 1; i < BLOCK_SIZE; i++) {
			tmp_block_write[k][i] = tmp_block_read[col_order[k]][i];
		}
		tmp_block_write[k][k] = tmp_block_read[col_order[k]][k];

		// For each column right of current diagonal element
		for (int j = k + 1; j < BLOCK_SIZE; j++) {
			// For each element below it
			#pragma unroll
			for (int i = k+1; i < BLOCK_SIZE; i++) {
				tmp_block_write[j][i] = tmp_block_read[col_order[j]][i] - tmp_scale_col[j] * tmp_block_read[col_order[k]][i];
			}
		}
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
left_blocks_c2(local const DATA_TYPE top_block[BLOCK_SIZE][BLOCK_SIZE],
				local const DATA_TYPE current_block_in[BLOCK_SIZE][BLOCK_SIZE],
				local DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE],
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
							tmp_block_read2[j][i] - tmp_scale_col[i]
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
top_blocks_c3(local const DATA_TYPE left_block[BLOCK_SIZE][BLOCK_SIZE],
			  local const DATA_TYPE current_block_in[BLOCK_SIZE][BLOCK_SIZE],
			  local DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE],
			  const uint ipvt[BLOCK_SIZE]) {
	DATA_TYPE tmp_block3[BLOCK_SIZE][BLOCK_SIZE];


	for (int j = 0; j < BLOCK_SIZE; j++) {
		#pragma unroll
		for (int i = 0; i <  BLOCK_SIZE; i++) {
			current_block_out[j][i] = current_block_in[j][i];
		}
	}

	// For each diagonal element in left block
	for (int k=0; k < BLOCK_SIZE; k++) {
#ifdef DEBUG
		printf("A(%d):\n", k);
		for (int j = 0; j < BLOCK_SIZE; j++) {
			for (int i = 0; i <  BLOCK_SIZE; i++) {
				printf("%f, ", left_block[j][i]);
			}
			for (int i = 0; i <  BLOCK_SIZE; i++) {
				printf("%f, ", current_block_out[j][i]);
			}
			printf("\n");
		}
		printf("\n\n");
#endif
		uint col_order[BLOCK_SIZE];
		for (int i=0; i < BLOCK_SIZE; i++) {
			col_order[i] = i;
		}
		col_order[k] = ipvt[k];
		col_order[ipvt[k]] = k;

		DATA_TYPE new_k_row[BLOCK_SIZE];
		#pragma unroll
		for (int i = 0; i < BLOCK_SIZE; i++) {
			new_k_row[i] = current_block_out[col_order[k]][i];
		}
		// For each column in current block
		for (int j = k; j < BLOCK_SIZE; j++) {
			DATA_TYPE multiply = 0.0;
			// For each element below it
			if (j > k) {
				multiply = left_block[j][k];
			}
			#pragma unroll
			for (int i = 0; i < BLOCK_SIZE; i++) {
				tmp_block3[j][i] = current_block_out[col_order[j]][i]
					- multiply * current_block_out[col_order[k]][i];
			}
		}
		for (int j = k; j < BLOCK_SIZE; j++) {
			#pragma unroll
			for (int i = 0; i <  BLOCK_SIZE; i++) {
				current_block_out[j][i] = tmp_block3[j][i];
			}
		}
	}
	#ifdef DEBUG
		printf("A(end):\n");
		for (int j = 0; j < BLOCK_SIZE; j++) {
			for (int i = 0; i <  BLOCK_SIZE; i++) {
				printf("%f, ", left_block[j][i]);
			}
			for (int i = 0; i <  BLOCK_SIZE; i++) {
				printf("%f, ", current_block_out[j][i]);
			}
			printf("\n");
		}
		printf("\n\n");
	#endif
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
inner_blocks_c4(local const DATA_TYPE left_block[BLOCK_SIZE][BLOCK_SIZE],
				local const DATA_TYPE top_block[BLOCK_SIZE][BLOCK_SIZE],
				local const DATA_TYPE current_block_in[BLOCK_SIZE][BLOCK_SIZE],
				local DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE]) {
	DATA_TYPE tmp_block4[BLOCK_SIZE][BLOCK_SIZE];
	copy_block(current_block_in, tmp_block4);
	// For each diagonal element in left block

	for (int k=0; k < BLOCK_SIZE; k++) {
		// For each column in top block
		#pragma unroll
		for (int j = 0; j < BLOCK_SIZE; j++) {
			// For each element below it in current block
			#pragma unroll
			for (int i = 0; i < BLOCK_SIZE; i++) {
				current_block_out[j][i] = tmp_block4[j][i] - left_block[j][k] * top_block[k][i];
			}
		}
		copy_block(current_block_out, tmp_block4);
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
void gefa(global DATA_TYPE* restrict a, global uint* restrict pvt,  uint a_size) {

	local DATA_TYPE top_block[BLOCK_SIZE][BLOCK_SIZE];
	local DATA_TYPE left_block[BLOCK_SIZE][BLOCK_SIZE];
	local DATA_TYPE diag_block[BLOCK_SIZE][BLOCK_SIZE];
	local DATA_TYPE top_block_out[BLOCK_SIZE][BLOCK_SIZE];
	local DATA_TYPE left_block_out[BLOCK_SIZE][BLOCK_SIZE];
	local DATA_TYPE diag_block_out[BLOCK_SIZE][BLOCK_SIZE];
	local DATA_TYPE current_block[BLOCK_SIZE][BLOCK_SIZE];
	local DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE];


	// For each diagonal block do the following
	#pragma disable_loop_pipelining
	for (int diagonal_block=0; diagonal_block < a_size; diagonal_block++) {

		// load next block for factorization
		load_block(diag_block, a, diagonal_block, diagonal_block, a_size);

		DATA_TYPE scale_factors[BLOCK_SIZE];
		uint ipvt[BLOCK_SIZE];

		// execute factorization of next block
		lu_factorization_c1(diag_block, diag_block_out, scale_factors,
													ipvt);

		#pragma unroll GLOBAL_MEM_UNROLL
		for (int i=0; i<BLOCK_SIZE; i++) {
			pvt[diagonal_block * BLOCK_SIZE + i] = diagonal_block * BLOCK_SIZE
																+ ipvt[i];
		}

		store_block(diag_block_out, a, diagonal_block, diagonal_block, a_size);

		// process all blocks columnwise
		for (int inner_x_block = diagonal_block + 1; inner_x_block < a_size;
			inner_x_block++) {
				// update top block
				load_block(top_block, a, inner_x_block, diagonal_block, a_size);
				top_blocks_c3(diag_block_out, top_block, top_block_out, ipvt);
				store_block(top_block_out, a, inner_x_block,
														diagonal_block, a_size);

				for (int inner_y_block = diagonal_block + 1;
									inner_y_block < a_size; inner_y_block++) {

						// update left block, if it was not already done
						if (inner_x_block == diagonal_block + 1) {
							load_block(left_block, a, diagonal_block,
														inner_y_block, a_size);
							left_blocks_c2(diag_block_out, left_block,
												left_block_out, scale_factors);
							store_block(left_block_out, a, diagonal_block,
														inner_y_block, a_size);
						} else {
							load_block(left_block_out, a, diagonal_block,
														inner_y_block, a_size);
						}

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
