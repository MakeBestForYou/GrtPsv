#include "scularfunction.h"
using namespace Eigen;
using namespace std;

ScularFunction::ScularFunction(model* media,char* scular_file)
{
    this->media = media;
    if(!(this->scular_out_fp = fopen(scular_file,"w+")))
    {
        printf("can't open scular function file!/n");
        exit(1);
    }
    this->E = this->malloc_array3(this->media->layer,4,4);
    this->RT = this->malloc_array3(this->media->layer,4,4);
    this->A_u = this->malloc_array3(this->media->layer,2,2);
    this->A_d = this->malloc_array3(this->media->layer,2,2);
    this->v = this->malloc_array1_cf(this->media->N_k);
    this->y = this->malloc_array1_cf(this->media->N_k);
}

void ScularFunction::calculate_scular(float frequency)
{
    for(int i = 0;i<this->media->N_k;i++)
    {
        this->calculate_E(frequency,this->media->k[i]);
    }
}

void ScularFunction::calculate_R_T(float frequency, float wavenumber)
{
    this->calculate_E(frequency,wavenumber);
    float** temp = this->malloc_array2(4,4);

    Matrix<float,2,2> E_11;
    Matrix<float,2,2> E_12;
    Matrix<float,2,2> E_21;
    Matrix<float,2,2> E_22;
    MatrixXf temp_matrix(4,4);
    MatrixXf A_matrix_0(4,4);
    MatrixXf A_matrix_1(4,4);
    MatrixXf E_matrix_0(4,4);
    MatrixXf E_matrix_1(4,4);
    MatrixXf E_matrix_2(4,4);
    MatrixXf E_matrix_3(4,4);
    MatrixXf RT_matrix(4,4);

    for(int j = 1;j < media->layer - 1;j++)
    {
        this->calculate_A(media->depth[j]);

        for(int ii = 0;ii<4;ii++)
            for(int jj = 0;jj<4;jj++)
            {
                E_matrix_0(ii,jj) = E[j][ii][jj];
                E_matrix_1(ii,jj) = E[j-1][ii][jj];

                A_matrix_0(ii,jj) = A[j][ii][jj];
                A_matrix_1(ii,jj) = A[j-1][ii][jj];
            }

        E_matrix_2.block(0,0,2,2) = E_matrix_0.block(0,0,2,2);
        E_matrix_2.block(0,2,2,2) = -E_matrix_1.block(0,2,2,2);
        E_matrix_2.block(2,0,2,2) = E_matrix_0.block(2,0,2,2);
        E_matrix_2.block(2,2,2,2) = -E_matrix_1.block(2,2,2,2);

        E_matrix_3.block(0,0,2,2) = E_matrix_1.block(0,0,2,2);
        E_matrix_3.block(0,2,2,2) = -E_matrix_0.block(0,2,2,2);
        E_matrix_3.block(2,0,2,2) = E_matrix_1.block(2,0,2,2);
        E_matrix_3.block(2,2,2,2) = -E_matrix_0.block(2,2,2,2);

        temp_matrix.block(0,2,2,2) = MatrixXf::Zero(2,2);
        temp_matrix.block(2,0,2,2) = MatrixXf::Zero(2,2);
        temp_matrix.block(0,0,2,2) = A_matrix_1.block(0,0,2,2);
        temp_matrix.block(2,2,2,2) = A_matrix_0.block(0,0,2,2);

        RT_matrix = E_matrix_0.inverse()*E_matrix_1*temp_matrix;

        for(int ii = 0;ii<4;ii++)
            for(int jj = 0;jj<4;jj++)
            {
                this->RT[j][ii][jj] = RT_matrix(ii,jj);
            }
    }
    int j = media->layer - 1;
    this->calculate_A(media->depth[j]);

    for(int ii = 0;ii<4;ii++)
        for(int jj = 0;jj<4;jj++)
        {
            E_matrix_0(ii,jj) = E[j][ii][jj];
            E_matrix_1(ii,jj) = E[j-1][ii][jj];

            A_matrix_0(ii,jj) = A[j][ii][jj];
            A_matrix_1(ii,jj) = A[j-1][ii][jj];
        }

    E_matrix_2.block(0,0,2,2) = E_matrix_0.block(0,0,2,2);
    E_matrix_2.block(0,2,2,2) = -E_matrix_1.block(0,2,2,2);
    E_matrix_2.block(2,0,2,2) = E_matrix_0.block(2,0,2,2);
    E_matrix_2.block(2,2,2,2) = -E_matrix_1.block(2,2,2,2);

    temp_matrix.block(0,0,2,2) = E_matrix_1.block(0,0,2,2)*A_matrix_1.block(0,0,2,2);
    temp_matrix.block(2,0,2,2) = E_matrix_1.block(2,0,2,2)*A_matrix_1.block(0,0,2,2);

    MatrixXf temp1(4,4);
    temp1.block(0,0,4,4) = MatrixXf::Zero(4,4);
    temp1.block(0,0,4,2) = E_matrix_2.inverse() * temp_matrix.block(0,0,4,2);

    for(int ii = 0;ii<4;ii++)
        for(int jj = 0;jj<4;jj++)
        {
            this->RT[j][ii][jj] = temp1(ii,jj);
        }
}

void ScularFunction::calculate_A(float depth_z)
{
    float z = depth_z;
    for(int i = 0;i<media->layer - 1;i++)
    {
        this->A_d[i][0][0] = exp(-this->y[i]*(z - media->depth[i]));
        this->A_d[i][0][1] = 0.0;
        this->A_d[i][1][0] = 0.0;
        this->A_d[i][1][1] = exp(-this->v[i]*(z - media->depth[i]));

        this->A_u[i][0][0] = exp(-this->y[i]*(media->depth[i+1] - z));
        this->A_u[i][0][1] = 0.0;
        this->A_u[i][1][0] = 0.0;
        this->A_u[i][1][1] = exp(-this->v[i]*(media->depth[i+1] - z));
    }
    int i = media->layer - 1;
    this->A_d[i][0][0] = exp(-this->y[i]*(z - media->depth[i]));
    this->A_d[i][1][1] = exp(-this->v[i]*(z - media->depth[i]));
    for(int layer=0;layer<4;layer++)
        for(int i=0;i<2;i++)
            for(int j = 0;j<2;j++)
        {
            this->A[layer][i][j] = this->A_d[layer][i][j];
            this->A[layer][i+2][j+2] = this->A_u[layer][i][j];
        }
}

void ScularFunction::calculate_v_gamma(float frequency, float wavenumber)
{
    float w = 2*M_PI*frequency;

    for(int i = 0;i<this->media->layer;i++)
    {
        double tmp = pow(wavenumber,2);
        complex<double> wave_complex(pow(wavenumber,2),0.0);
        wave_complex = pow(wave_complex*4.0,0.5);
        this->v[i] = pow(pow(wavenumber,2) - pow(w/this->media->belta[i],2),0.5);
        this->y[i] = pow(pow(wavenumber,2) - pow(w/this->media->alpha[i],2),0.5);
    }
}

void ScularFunction::calculate_E(float frequency,float wavenumber)
{
    this->calculate_v_gamma(frequency,wavenumber);
    float w = 2*M_PI*frequency;
    float* X = this->malloc_array1(this->media->layer);
    for(int i = 0;i<this->media->layer;i++)
    {
        X[i] = wavenumber*wavenumber + pow(this->v[i],2);
    }

    for(int i = 0;i<this->media->layer;i++)
    {
        this->E[i][0][0] = media->alpha[i] * wavenumber / w;
        this->E[i][0][1] = media->belta[i] * this->v[i] / w;
        this->E[i][0][2] = media->alpha[i] * wavenumber / w;
        this->E[i][0][3] = media->belta[i] * this->v[i] / w;

        this->E[i][1][0] = media->alpha[i] * this->y[i] / w;
        this->E[i][1][1] = media->belta[i] * wavenumber / w;
        this->E[i][1][2] = -media->alpha[i] * this->y[i] / w;
        this->E[i][1][3] = -media->belta[i] * wavenumber / w;

        this->E[i][2][0] = -2*media->alpha[i] * media->miu[i] * wavenumber * this->y[i] / w;
        this->E[i][2][1] = -media->belta[i] * media->miu[i] * X[i] / w;
        this->E[i][2][2] = 2*media->alpha[i] * media->miu[i] * wavenumber * this->y[i] / w;
        this->E[i][2][3] = media->belta[i] * media->miu[i] * X[i] / w;

        this->E[i][3][0] = - media->alpha[i] * media->miu[i] * X[i] / w;
        this->E[i][3][1] = -2*media->belta[i] * media->miu[i] * wavenumber * this->v[i] / w;
        this->E[i][3][2] = -media->alpha[i] * media->miu[i] * X[i] / w;
        this->E[i][3][3] = -2*media->belta[i] * media->miu[i] * wavenumber * this->v[i] / w;
    }

    this->release_array1(X);
}

void ScularFunction::scular_output()
{
    for(int i = 0;i<this->scular_N;i++)
    {
        fprintf(this->scular_out_fp,"%f/t",this->scular[i]);
    }
    fprintf(this->scular_out_fp,"/n");
}

float*** ScularFunction::malloc_array3(int page, int row, int low)
{
    float*** array3 = (float***)malloc(sizeof(float**)*page);
    for(int i = 0;i<page;i++)
    {
        array3[i] = (float**)malloc(sizeof(float*)*row);
        for(int ii = 0;ii<4;ii++)
        {
            array3[i][ii] = (float*)malloc(sizeof(float)*low);
        }
    }
    return array3;
}



float** ScularFunction::malloc_array2(int row, int low)
{
    float** array2 = (float**)malloc(sizeof(float*)*row);
    for(int i = 0;i<row;i++)
    {
        array2[i] = (float*)malloc(sizeof(float)*low);
    }
    return array2;
}

float* ScularFunction::malloc_array1(int low)
{
    float* array1 = (float*)malloc(sizeof(float)*low);
    return array1;
}

complex<float>*** ScularFunction::malloc_array3_cf(int page, int row, int low)
{
    complex<float>*** array3 = (complex<float>***)malloc(sizeof(complex<float>**)*page);
    for(int i = 0;i<page;i++)
    {
        array3[i] = (complex<float>**)malloc(sizeof(complex<float>*)*row);
        for(int ii = 0;ii<4;ii++)
        {
            array3[i][ii] = (complex<float>*)malloc(sizeof(complex<float>)*low);
        }
    }
    return array3;
}

complex<float>** ScularFunction::malloc_array2_cf(int row, int low)
{
    complex<float>** array2 = (complex<float>**)malloc(sizeof(complex<float>*)*row);
    for(int i = 0;i<row;i++)
    {
        array2[i] = (complex<float>*)malloc(sizeof(complex<float>)*low);
    }
    return array2;
}

complex<float>* ScularFunction::malloc_array1_cf(int low)
{
    complex<float>* array1 = (complex<float>*)malloc(sizeof(complex<float>)*low);
    return array1;
}

void ScularFunction::release_array3(float ***array3, int page, int row)
{
    for(int i = 0;i<page;i++)
    {
        for(int ii = 0;ii<row;ii++)
        {
            free(array3[i][ii]);
        }
    }
}

void ScularFunction::release_array2(float **array2, int row)
{
    for(int i = 0;i<row;i++)
    {
        free(array2[i]);
    }
}

void ScularFunction::release_array1(float *array1)
{
    free(array1);
}

void ScularFunction::release_array3_cf(complex<float> ***array3, int page, int row)
{
    for(int i = 0;i<page;i++)
    {
        for(int ii = 0;ii<row;ii++)
        {
            free(array3[i][ii]);
        }
    }
}

void ScularFunction::release_array2_cf(complex<float>** array2, int row)
{
    for(int i = 0;i<row;i++)
    {
        free(array2[i]);
    }
}

void ScularFunction::release_array1_cf(complex<float>* array1)
{
    free(array1);
}


