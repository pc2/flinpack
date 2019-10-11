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
Copy a block to another memory location column-wise

@param a_from local memory buffer to load the block from
@param a_to the local memory buffer to copy the block into
*/
void
copy_block_col(local const DATA_TYPE a_from[BLOCK_SIZE][BLOCK_SIZE],
				DATA_TYPE a_to[BLOCK_SIZE][BLOCK_SIZE]) {

	for (int i = 0; i < BLOCK_SIZE; i++) {
		#pragma unroll
		for (int j = 0; j < BLOCK_SIZE; j++) {
			a_to[j][i] = a_from[j][i];
		}
	}
}


/**
Standard LU factorization on a block with fixed size

Case 1 of Zhangs description
*/
void
lu_factorization_c1(local const DATA_TYPE a_block_in[BLOCK_SIZE][BLOCK_SIZE],
					local DATA_TYPE a_block_out[BLOCK_SIZE][BLOCK_SIZE],
					DATA_TYPE scale_factors[BLOCK_SIZE]) {
	// copy columnwise
	for (int i = 0; i < BLOCK_SIZE; i++) {
		for (int j = 0; j <  BLOCK_SIZE; j++) {
			a_block_out[j][i] = a_block_in[j][i];
		}
	}
	// For each diagnonal element
	for (int k = 0; k < BLOCK_SIZE; k++) {
		DATA_TYPE tmp_scale_col[BLOCK_SIZE];
		scale_factors[k] = 1.0 / a_block_out[k][k];
		// For each element below it
		for (int i = k + 1; i < BLOCK_SIZE; i++) {
			tmp_scale_col[i] =  a_block_out[i][k] * scale_factors[k];
			a_block_out[i][k] = tmp_scale_col[i];
		}
		// For each column right of current diagonal element
		for (int j = k + 1; j < BLOCK_SIZE; j++) {
			// For each element below it
			#pragma unroll
			for (int i = k+1; i < BLOCK_SIZE; i++) {
				DATA_TYPE tmp = a_block_out[i][j] - tmp_scale_col[i] * a_block_out[k][j];
				a_block_out[i][j] = tmp;
			}
		}
	}
}

/**
Modifying the blocks on the leftmost side

Case 2 of Zhangs description
*/
void
left_blocks_c2(local const DATA_TYPE top_block[BLOCK_SIZE][BLOCK_SIZE],
				local const DATA_TYPE current_block_in[BLOCK_SIZE][BLOCK_SIZE],
				local DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE],
				const DATA_TYPE scale_factors[BLOCK_SIZE]) {
	DATA_TYPE tmp_block2[BLOCK_SIZE][BLOCK_SIZE];
	DATA_TYPE tmp_scale_col[BLOCK_SIZE];
	copy_block_col(current_block_in, tmp_block2);
	// For each diagonal element in top block
	for (int k=0; k < BLOCK_SIZE; k++) {
		// For each element below it in current block
		#pragma unroll
		for (int i=0; i < BLOCK_SIZE; i++) {
			// printf("C2: %f * %f\n",tmp_block2[i][k], scale_factors[k]);
			tmp_scale_col[i] = tmp_block2[i][k] * scale_factors[k];
			current_block_out[i][k] = tmp_scale_col[i];
			tmp_block2[i][k] = tmp_scale_col[i];
		}
		// For each column right of the current diagnonal element
		#pragma max_concurrency 1
		for (int i = k+1; i < BLOCK_SIZE; i++) {
			DATA_TYPE tmp_col[BLOCK_SIZE];
			#pragma unroll
			for (int j = 0; j < BLOCK_SIZE; j++) {
				tmp_col[j] = tmp_block2[j][i] - tmp_scale_col[j] * top_block[k][i];
				current_block_out[j][i] = tmp_col[j];
			}
			#pragma unroll
			for (int j = 0; j < BLOCK_SIZE; j++) {
				tmp_block2[j][i] = tmp_col[j];
			}
		}
	}
}

/**
Modifying the blocks on the top but not on the left

Case 3 of Zhangs description
*/
void
top_blocks_c3(local const DATA_TYPE left_block[BLOCK_SIZE][BLOCK_SIZE],
			  local const DATA_TYPE current_block_in[BLOCK_SIZE][BLOCK_SIZE],
			  local DATA_TYPE current_block_out[BLOCK_SIZE][BLOCK_SIZE]) {
	DATA_TYPE tmp_block3[BLOCK_SIZE][BLOCK_SIZE];
	copy_block(current_block_in, tmp_block3);
	for (int i = 0; i < BLOCK_SIZE; i++) {
		current_block_out[0][i] = current_block_in[0][i];
	}
	// For each diagonal element in left block
	#pragma loop_coalesce 2
	for (int k=0; k < BLOCK_SIZE; k++) {
		// For each column in current block
		for (int j = k+1; j < BLOCK_SIZE; j++) {
			// For each element below it
			DATA_TYPE tmp_row[BLOCK_SIZE];
			#pragma unroll
			for (int i = 0; i < BLOCK_SIZE; i++) {
				tmp_row[i] = tmp_block3[j][i] - left_block[j][k] * tmp_block3[k][i];
				current_block_out[j][i] = tmp_row[i];
			}
			#pragma unroll
			for (int i = 0; i < BLOCK_SIZE; i++) {
				tmp_block3[j][i] = tmp_row[i];
			}
		}
	}
}

/**
Modifying the inner blocks

Case 4 of Zhangs description
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
@param a_size the x and y size of the matrix in blocks
*/
__attribute__((uses_global_work_offset(0)))
__kernel
void gefa(global DATA_TYPE* restrict a, uint a_size) {

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

		// execute factorization of next block
		lu_factorization_c1(diag_block, diag_block_out, scale_factors);

		store_block(diag_block_out, a, diagonal_block, diagonal_block, a_size);

		for (int inner_x_block = diagonal_block + 1; inner_x_block < a_size;
			inner_x_block++) {
				// update top block
				load_block(top_block, a, inner_x_block, diagonal_block, a_size);
				top_blocks_c3(diag_block_out, top_block, top_block_out);
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
