#include <stdio.h>
#include <stdlib.h>
#include <omp.h>  // For future parallel version
#include <time.h>

// Function to read Matrix Market file and load it into COO format
void read_matrix_market(const char *filename, int *nrows, int *ncols, int *nnz, 
                        int **row_idx, int **col_idx, double **values) {
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        printf("Error: Unable to open file %s\n", filename);
        exit(1);
    }

    // Skip comments and read header
    char line[256];
    do {
        fgets(line, sizeof(line), f);
    } while (line[0] == '%');

    // Read dimensions and number of non-zeros
    sscanf(line, "%d %d %d", nrows, ncols, nnz);

    // Allocate memory for COO representation
    *row_idx = (int *)malloc((*nnz) * sizeof(int));
    *col_idx = (int *)malloc((*nnz) * sizeof(int));
    *values  = (double *)malloc((*nnz) * sizeof(double));

    // Read non-zero entries (COO format)
    for (int i = 0; i < *nnz; i++) {
        fscanf(f, "%d %d %lf", &(*row_idx)[i], &(*col_idx)[i], &(*values)[i]);
        // Convert 1-based to 0-based indexing
        (*row_idx)[i]--;
        (*col_idx)[i]--;
    }

    fclose(f);
}

// Function to convert COO to CSR format
void convert_to_csr(int nrows, int ncols, int nnz, 
                    int *row_idx, int *col_idx, double *values, 
                    int **row_ptr, int **csr_col_idx, double **csr_values) {
    // Allocate memory for CSR arrays
    *row_ptr = (int *)malloc((nrows + 1) * sizeof(int));
    *csr_col_idx = (int *)malloc(nnz * sizeof(int));
    *csr_values = (double *)malloc(nnz * sizeof(double));

    // Initialize row_ptr with 0s
    for (int i = 0; i <= nrows; i++) {
        (*row_ptr)[i] = 0;
    }

    // Count the number of non-zeros per row
    for (int i = 0; i < nnz; i++) {
        (*row_ptr)[row_idx[i] + 1]++;
    }

    // Cumulative sum to get the starting index for each row
    for (int i = 0; i < nrows; i++) {
        (*row_ptr)[i + 1] += (*row_ptr)[i];
    }

    // Fill in the csr_col_idx and csr_values arrays
    int *current_row_position = (int *)malloc(nrows * sizeof(int));
    for (int i = 0; i < nrows; i++) {
        current_row_position[i] = (*row_ptr)[i];
    }

    for (int i = 0; i < nnz; i++) {
        int row = row_idx[i];
        int idx = current_row_position[row]++;
        (*csr_col_idx)[idx] = col_idx[i];
        (*csr_values)[idx] = values[i];
    }

    free(current_row_position);
}

// Define a tile size for CSR5 (e.g., 32 rows per tile)
#define TILE_SIZE 32

// Structure to hold CSR5 metadata for each tile
typedef struct {
    int tile_start;      // Starting index of the tile in csr_values
    int num_nonzeros;    // Total number of non-zeros in the tile
} CSR5Tile;

// Function to convert CSR to CSR5 format
void convert_to_csr5(int nrows, int *row_ptr, int nnz, int *csr_col_idx, double *csr_values, 
                     CSR5Tile **csr5_tiles, int *num_tiles) {
    *num_tiles = (nrows + TILE_SIZE - 1) / TILE_SIZE;  // Calculate number of tiles
    *csr5_tiles = (CSR5Tile *)malloc((*num_tiles) * sizeof(CSR5Tile));

    int tile_idx = 0;
    for (int row = 0; row < nrows; row += TILE_SIZE) {
        int tile_start = row_ptr[row];  // First non-zero element in this tile
        int tile_end = (row + TILE_SIZE < nrows) ? row_ptr[row + TILE_SIZE] : nnz;  // End of the tile
        int num_nonzeros = tile_end - tile_start;

        (*csr5_tiles)[tile_idx].tile_start = tile_start;
        (*csr5_tiles)[tile_idx].num_nonzeros = num_nonzeros;

        tile_idx++;
    }
}

// Function to perform SpMV using CSR5 format
void spmv_csr5(int nrows, CSR5Tile *csr5_tiles, int *row_ptr, int *csr_col_idx, double *csr_values, 
               double *vector, double *result, int num_tiles) {
    for (int tile = 0; tile < num_tiles; tile++) {
        int tile_start = csr5_tiles[tile].tile_start;
        int num_nonzeros = csr5_tiles[tile].num_nonzeros;

        // Process each row in the current tile
        for (int row = tile * TILE_SIZE; row < (tile + 1) * TILE_SIZE && row < nrows; row++) {
            double sum = 0.0;
            for (int idx = row_ptr[row]; idx < row_ptr[row + 1]; idx++) {
                sum += csr_values[idx] * vector[csr_col_idx[idx]];
            }
            result[row] = sum;  // Store the result in the output vector
        }
    }
}

void read_vector_market(const char *filename, int nrows, double *vector) {
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        printf("Error: Unable to open vector file %s\n", filename);
        exit(1);
    }

    // Skip comments and read header
    char line[256];
    do {
        fgets(line, sizeof(line), f);
    } while (line[0] == '%');

    // Read vector entries (MM format assumes 1 column for vectors)
    for (int i = 0; i < nrows; i++) {
        int row;
        double value;
        fscanf(f, "%d %lf", &row, &value);
        // Convert 1-based to 0-based indexing
        vector[row - 1] = value;
    }

    fclose(f);
}

void write_vector_market(const char *filename, int nrows, double *vector) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error: Unable to open output file %s\n", filename);
        exit(1);
    }

    // Write the header
    fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%d 1\n", nrows);

    // Write the vector values
    for (int i = 0; i < nrows; i++) {
        fprintf(f, "%d %.6f\n", i + 1, vector[i]);
    }

    fclose(f);
}


int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s matrix.mtx vector.mtx output.mtx\n", argv[0]);
        return 1;
    }

    const char *matrix_file = argv[1];
    const char *vector_file = argv[2];
    const char *output_file = argv[3];

    int nrows, ncols, nnz;
    int *row_idx, *col_idx;
    double *values;

    // Step 1: Read the matrix in COO format
    read_matrix_market(matrix_file, &nrows, &ncols, &nnz, &row_idx, &col_idx, &values);

    // Step 2: Convert to CSR format
    int *row_ptr, *csr_col_idx;
    double *csr_values;
    convert_to_csr(nrows, ncols, nnz, row_idx, col_idx, values, &row_ptr, &csr_col_idx, &csr_values);

    // Step 3: Convert CSR to CSR5 format
    CSR5Tile *csr5_tiles;
    int num_tiles;
    convert_to_csr5(nrows, row_ptr, nnz, csr_col_idx, csr_values, &csr5_tiles, &num_tiles);

    // Step 4: Load the vector from vector.mtx (Assume it is a single column vector)
    double *vector = (double *)malloc(ncols * sizeof(double));  // Input vector
    double *result = (double *)malloc(nrows * sizeof(double));  // Output vector
    read_vector_market(vector_file, ncols, vector);

    // Start timing
    clock_t start_time = clock();

    // Step 5: Perform SpMV using CSR5 format
    spmv_csr5(nrows, csr5_tiles, row_ptr, csr_col_idx, csr_values, vector, result, num_tiles);

    // End timing
    clock_t end_time = clock();

    // Calculate elapsed time in seconds
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("Execution Time: %.6f seconds\n", elapsed_time);

    // Step 6: Write the result vector to output.mtx
    write_vector_market(output_file, nrows, result);

    // Clean up memory
    free(row_idx);
    free(col_idx);
    free(values);
    free(row_ptr);gcc -o seq_spmv_csr5 seq_spmv_csr5.c

    free(csr_col_idx);
    free(csr_values);
    free(csr5_tiles);
    free(vector);
    free(result);

    return 0;
}
