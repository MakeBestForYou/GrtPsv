#ifndef SCULARFUNCTION_H
#define SCULARFUNCTION_H
#include"model.h"
#include<stdlib.h>
#include<math.h>
#include<complex>
#include<eigen3/Eigen/Dense>
using namespace std;

class ScularFunction
{
public:
    ScularFunction(model* media,char* output_file);
    void calculate_scular(float frequency);
    void scular_output();

public:
    float* scular;
    float* velocity;
    int scular_N;

private:
    model* media;
    FILE* scular_out_fp;
    float*** E;

    complex<float>* v;
    complex<float>* y;

    float*** RT;
    float*** A;
    float*** A_u;
    float*** A_d;

    float* R_ud;
    float* R_du;
    float* T_d;
    float* T_u;

private:
    float*** malloc_array3(int page,int row,int low);
    float** malloc_array2(int row,int low);
    float* malloc_array1(int low);

    void release_array3(float*** array3,int page,int row);
    void release_array2(float** array2,int row);
    void release_array1(float* array1);

    complex<float>*** malloc_array3_cf(int page,int row,int low);
    complex<float>** malloc_array2_cf(int row,int low);
    complex<float>* malloc_array1_cf(int low);

    void release_array3_cf(complex<float>*** array3,int page,int row);
    void release_array2_cf(complex<float>** array2,int row);
    void release_array1_cf(complex<float>* array1);

    void calculate_R_T(float frequency,float wavenumber);
    void calculate_E(float frequency,float wavenumber);
    void calculate_v_gamma(float frequency,float wavenumber);
    void calculate_A(float depth_z);

};

#endif // SCULARFUNCTION_H
