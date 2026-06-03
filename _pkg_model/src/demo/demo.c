/*
 *  This is a set of testbench code for handling large dynamically allocated arrays in C
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "complex.hpp"

void accept_vector(double *x) {
    /* Print simple confirmation that everything is working */
    printf("Yo");
}   

void accept_matrix(double **x) {
    /* Print simple confirmation that everything is working */
    printf("Yo");
    printf("\n");

    /* Loop through elements and print */
    for (int i=0; i<2; i++) {
        for (int j=0; j<3; j++) {
            printf("%.6f", x[i][j]);
            printf("\n");
        }
    }
}   

void allocate_vector_complex() {
    /* Allocate memory for a complex vector */
    COMPLEX *x;
    x = (COMPLEX*) calloc(4, sizeof(COMPLEX));

    for (int i = 0; i < 4; i++) {
        printf("x = %f\n", x[i].x);
        printf("y = %f\n", x[i].y);
    }

    /* Free memory */
    free(x);
}

void allocate_matrix_complex() {
    /* Declare print_vector */
    void print_vector(COMPLEX *x);

    /* Allocate memory for a complex vector */
    COMPLEX **x;
    for (int i = 0; i < 5; i++) {
        x[i] = (COMPLEX*) calloc(4, sizeof(COMPLEX));
    }

    // for (int i = 0; i < 5; i++) {
    //     for (int j = 0; j < 4; j++){
    //         printf("x = %f\n", x[i][j].x);
    //         printf("y = %f\n", x[i][j].y);
    //     }
    // }

    print_vector(x[1]);

    /* Free memory */
    for (int i = 0; i < 5; i++) {
        free(x[i]);
    }
}

void allocate_matrix_complex2() {
    /* Declare print_vector */
    void print_vector(COMPLEX *x);

    /* Allocate memory for a complex vector */
    COMPLEX *x[5];
    for (int i = 0; i < 5; i++) {
        x[i] = (COMPLEX*) calloc(4, sizeof(COMPLEX));
    }

    // for (int i = 0; i < 5; i++) {
    //     for (int j = 0; j < 4; j++){
    //         printf("x = %f\n", x[i][j].x);
    //         printf("y = %f\n", x[i][j].y);
    //     }
    // }

    print_vector(x[1]);

    /* Free memory */
    for (int i = 0; i < 5; i++) {
        free(x[i]);
    }
}


void print_vector(COMPLEX *x) {
    for (int i = 0; i < 4; i++) {
            printf("x = %f\n", x[i].x);
            printf("y = %f\n", x[i].y);
    }
}

void test_wbgt(int *x) {
    printf("Original value = %d\n", *x);
    (*x)++;
}

void pass_by_reference() {
    int vector_of_ints[10];
    vector_of_ints[3] = 22;
    test_wbgt(&vector_of_ints[3]);
    printf("New value = %d\n", vector_of_ints[3]);
}

void allocate_3d_array() {
    int M = 10;
    int N = 4;
    int O = 12;

    int*** test = (int***) malloc(M * sizeof(int**));

    for (int i = 0; i < M; i++) {
        test[i] = (int**) malloc(N * sizeof(int*));
        for (int j = 0; j < N; j++) {
            test[i][j] = (int*) malloc(O * sizeof(int));
        }
    }

    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 12; k++) {
                test[i][j][k] = 10*i + j - k;
            }
        }
    }
    printf("Val: %d\n", test[1][2][3]);
}

void modify_vector_of_matrix_inplace() {
    double* tmp[10];
    for (int i = 0; i < 10; i++) {
        tmp[i] = (double*) calloc(500, sizeof(double));
    }
    modify(tmp[3]);
    printf("Modified value is %f\n", tmp[3][3]);
}

void modify(double *x) {
    x[3] = 3.0;
}