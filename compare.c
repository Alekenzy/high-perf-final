#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define precision for comparison
#define FLOAT_TOLERANCE 1e-5   // For single precision comparison

// Function to read a vector from a Matrix Market format file
void read_vector_market(const char *filename, int *nrows, double **vector) {
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

    // Read the number of rows (the vector has 1 column)
    fscanf(f, "%d %*d", nrows);

    // Allocate memory for the vector
    *vector = (double *)malloc((*nrows) * sizeof(double));

    // Read vector entries
    for (int i = 0; i < *nrows; i++) {
        int row;
        double value;
        fscanf(f, "%d %lf", &row, &value);
        (*vector)[row - 1] = value;  // 1-based to 0-based index
    }

    fclose(f);
}

// Function to compare two vectors with a tolerance
int compare_vectors(int nrows1, double *vector1, int nrows2, double *vector2, double tolerance) {
    if (nrows1 != nrows2) {
        printf("Vector lengths differ: %d vs %d\n", nrows1, nrows2);
        return 0;  // Fail
    }

    // Compare each element
    for (int i = 0; i < nrows1; i++) {
        if (fabs(vector1[i] - vector2[i]) > tolerance) {
            printf("Difference at index %d: %f vs %f\n", i, vector1[i], vector2[i]);
            return 0;  // Fail
        }
    }

    return 1;  // Pass
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s output1.mtx output2.mtx\n", argv[0]);
        return 1;
    }

    const char *output1_file = argv[1];
    const char *output2_file = argv[2];

    int nrows1, nrows2;
    double *vector1, *vector2;

    // Step 1: Read both result vectors
    read_vector_market(output1_file, &nrows1, &vector1);
    read_vector_market(output2_file, &nrows2, &vector2);

    // Step 2: Compare the vectors
    int result = compare_vectors(nrows1, vector1, nrows2, vector2, FLOAT_TOLERANCE);

    // Step 3: Output "Pass" or "Fail"
    if (result) {
        printf("Pass\n");
    } else {
        printf("Fail\n");
    }

    // Clean up memory
    free(vector1);
    free(vector2);

    return 0;
}
