#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void generate_vector_market(const char *filename, int nrows) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error: Unable to open file %s\n", filename);
        exit(1);
    }

    // Write the header for Matrix Market format
    fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%d 1\n", nrows);

    // Write the vector values
    for (int i = 1; i <= nrows; i++) {
        fprintf(f, "%d %.6f\n", i, (double)rand() / RAND_MAX);
    }

    fclose(f);
}

int get_matrix_dimensions(const char *filename, int *nrows, int *ncols) {
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        printf("Error: Unable to open file %s\n", filename);
        return 1;
    }

    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '%') {
            continue;  // Skip comment lines
        }

        // Read matrix dimensions from the first non-comment line
        if (sscanf(line, "%d %d", nrows, ncols) == 2) {
            fclose(f);
            return 0;  // Successfully read dimensions
        }
    }

    fclose(f);
    printf("Error: Unable to read matrix dimensions from %s\n", filename);
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc < 1) {
        printf("Usage: %s matrix.mtx \n", argv[0]);
        return 1;
    }
    
    const char *matrix_filename = argv[1];
    int nrows, ncols;

    if (get_matrix_dimensions(matrix_filename, &nrows, &ncols) == 0) {
        printf("Matrix dimensions: %d x %d\n", nrows, ncols);
        // Generate vector file of the same size as the matrix's number of rows
        char vector_filename[256];
        snprintf(vector_filename, sizeof(vector_filename), "vector_%s", matrix_filename);
        generate_vector_market(vector_filename, nrows);  // Generate vector of size nrows
        printf("Vector file generated: %s\n", vector_filename);
    } else {
        printf("Failed to retrieve matrix dimensions.\n");
    }

    return 0;
}
