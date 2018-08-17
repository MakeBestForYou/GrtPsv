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
    this->scular = this->malloc_array1_cf(this->media->N_k);
    this->velocity = this->malloc_array1(this->media->N_k);
    this->E = this->malloc_array3_cf(this->media->layer,4,4);
    this->RT = this->malloc_array3_cf(this->media->layer,4,4);
    this->R_du_g = this->malloc_array3_cf(this->media->layer,2,2);
    this->R_ud_g = this->malloc_array3_cf(this->media->layer,2,2);
    this->A = this->malloc_array3_cf(this->media->layer,4,4);
    this->A_u = this->malloc_array3_cf(this->media->layer,2,2);
    this->A_d = this->malloc_array3_cf(this->media->layer,2,2);
    this->v = this->malloc_array1_cf(this->media->layer);
    this->y = this->malloc_array1_cf(this->media->layer);
}

void ScularFunction::realef_array()
{
    this->release_array1_cf(this->scular);
    this->release_array1(this->velocity);
    this->release_array3_cf(this->E,this->media->layer,4);
    this->release_array3_cf(this->RT,this->media->layer,4);
    this->release_array3_cf(this->R_du_g,this->media->layer,2);
    this->release_array3_cf(this->R_ud_g,this->media->layer,2);
    this->release_array3_cf(this->A,this->media->layer,4);
    this->release_array3_cf(this->A_u,this->media->layer,2);
    this->release_array3_cf(this->A_d,this->media->layer,2);
    this->release_array1_cf(this->v);
    this->release_array1_cf(this->y);
}

void ScularFunction::calculate_scular(float frequency)
{
    for(int i = 0;i<this->media->N_k;i++)
    {
        this->calculate_R_T(frequency,this->media->k[i]);
        this->calculate_RT_g();
        this->calculate_Det(i);
    }
}

complex<float> ScularFunction::calculate_scular(float frequency,float wavenumber)
{
    this->calculate_R_T(frequency,wavenumber);
    this->calculate_RT_g();
    return this->calculate_Det();
}

complex<float> ScularFunction::calculate_Det()
{
    MatrixXcf R_ud_g_(2,2);
    MatrixXcf R_du_g_(2,2);
    MatrixXcf dialog(2,2);
    dialog<<1,0,0,1;
    for(int ii = 0;ii<2;ii++)
        for(int jj = 0;jj<2;jj++)
        {
            R_ud_g_(ii,jj) = this->R_ud_g[0][ii][jj];
            R_du_g_(ii,jj) = this->R_du_g[1][ii][jj];
        }
    MatrixXcf tmp(2,2);
    tmp = dialog-R_ud_g_*R_du_g_;
    return tmp.determinant();
}

void ScularFunction::calculate_Det(int N_k)
{
    MatrixXcf R_ud_g_(2,2);
    MatrixXcf R_du_g_(2,2);
    MatrixXcf dialog(2,2);
    dialog<<1,0,0,1;
    for(int ii = 0;ii<2;ii++)
        for(int jj = 0;jj<2;jj++)
        {
            R_ud_g_(ii,jj) = this->R_ud_g[0][ii][jj];
            R_du_g_(ii,jj) = this->R_du_g[1][ii][jj];
        }
    MatrixXcf tmp(2,2);
    tmp = dialog-R_ud_g_*R_du_g_;
//    cout<<tmp.determinant()<<endl;
//    cout<<R_du_g<<R_ud_g<<endl;
    this->scular[N_k] =tmp.determinant();
}

void ScularFunction::calculate_RT_g()
{
    //first to compute the surface's genreal r-t cofficient
    this->calculate_A(media->depth[0]);
    MatrixXcf temp_1(2,2);
    MatrixXcf temp_2(2,2);
    MatrixXcf temp_3(2,2);
    MatrixXcf temp_4(2,2);
    MatrixXcf temp_5(2,2);
    MatrixXcf dialog(2,2);
    MatrixXcf result_R_du_g(2,2);
    dialog<<1,0,0,1;

    for(int ii = 0;ii<2;ii++)
        for(int jj = 0;jj<2;jj++)
        {
            temp_1(ii,jj) = E[1][ii+2][jj];
            temp_2(ii,jj) = E[1][ii+2][jj+2];
            temp_3(ii,jj) = A[1][ii+2][jj+2];
        }
    temp_4 = -temp_1.inverse()*temp_2*temp_3;
//    cout<<"R_ud_g"<<endl;
//    cout<<temp_4<<endl;
    for(int ii = 0;ii<2;ii++)
        for(int jj = 0;jj<2;jj++)
        {
            R_ud_g[0][ii][jj] = temp_4(ii,jj);
        }


    for(int i=0;i<media->layer;i++)
        for(int ii = 0;ii<2;ii++)
           for(int jj = 0;jj<2;jj++)
           {
                R_du_g[i][ii][jj] = RT[i][ii+2][jj];
            }

    for(int kk = media->layer-2;kk>=1;kk--)
    {
        for(int ii = 0;ii<2;ii++)
            for(int jj = 0;jj<2;jj++)
            {
                temp_1(ii,jj) = RT[kk][ii+2][jj];
                temp_2(ii,jj) = RT[kk][ii+2][jj+2];
                temp_4(ii,jj) = RT[kk][ii][jj+2];
                temp_5(ii,jj) = RT[kk][ii][jj];
                temp_3(ii,jj) = R_du_g[kk+1][ii][jj];
            }
        temp_4 = dialog-temp_4*temp_3;
        result_R_du_g = temp_1 + temp_2*temp_3*temp_4.inverse()*temp_5;
        for(int ii = 0;ii<2;ii++)
            for(int jj = 0;jj<2;jj++)
            {
                R_du_g[kk][ii][jj] = result_R_du_g(ii,jj);
            }
//        cout<<"R_ud_g"<<endl;
 //       cout<<result_R_du_g<<endl;
    }
}

void ScularFunction::calculate_R_T(float frequency, float wavenumber)
{
    this->calculate_E(frequency,wavenumber);
//    float** temp = this->malloc_array2(4,4);

/*    Matrix<float,2,2> E_11;
    Matrix<float,2,2> E_12;
    Matrix<float,2,2> E_21;
    Matrix<float,2,2> E_22;*/
    MatrixXcf temp_matrix(4,4);
    MatrixXcf A_matrix_0(4,4);
    MatrixXcf A_matrix_1(4,4);
    MatrixXcf E_matrix_0(4,4);
    MatrixXcf E_matrix_1(4,4);
    MatrixXcf E_matrix_2(4,4);
    MatrixXcf E_matrix_3(4,4);
    MatrixXcf RT_matrix(4,4);

    for(int j = 1;j < media->layer - 2;j++)
    {
        this->calculate_A(media->depth[j]);

        for(int ii = 0;ii<4;ii++)
            for(int jj = 0;jj<4;jj++)
            {
                E_matrix_0(ii,jj) = E[j][ii][jj];
                E_matrix_1(ii,jj) = E[j+1][ii][jj];

                A_matrix_0(ii,jj) = A[j][ii][jj];
                A_matrix_1(ii,jj) = A[j+1][ii][jj];
            }

        E_matrix_2.block(0,0,2,2) = E_matrix_0.block(0,0,2,2);
        E_matrix_2.block(0,2,2,2) = -E_matrix_1.block(0,2,2,2);
        E_matrix_2.block(2,0,2,2) = E_matrix_0.block(2,0,2,2);
        E_matrix_2.block(2,2,2,2) = -E_matrix_1.block(2,2,2,2);

        E_matrix_3.block(0,0,2,2) = E_matrix_1.block(0,0,2,2);
        E_matrix_3.block(0,2,2,2) = -E_matrix_0.block(0,2,2,2);
        E_matrix_3.block(2,0,2,2) = E_matrix_1.block(2,0,2,2);
        E_matrix_3.block(2,2,2,2) = -E_matrix_0.block(2,2,2,2);

        temp_matrix.block(0,2,2,2) = MatrixXcf::Zero(2,2);
        temp_matrix.block(2,0,2,2) = MatrixXcf::Zero(2,2);
        temp_matrix.block(0,0,2,2) = A_matrix_0.block(0,0,2,2);
        temp_matrix.block(2,2,2,2) = A_matrix_1.block(2,2,2,2);

//        cout<<temp_matrix<<endl;
//        cout<<E_matrix_2<<endl;
//        cout<<E_matrix_3<<endl;

        RT_matrix = E_matrix_3.inverse()*E_matrix_2*temp_matrix;
//        cout<<RT_matrix<<endl;
        for(int ii = 0;ii<4;ii++)
            for(int jj = 0;jj<4;jj++)
            {
                this->RT[j][ii][jj] = RT_matrix(ii,jj);
            }
    }
    int j = media->layer - 2;
    this->calculate_A(media->depth[j]);

    for(int ii = 0;ii<4;ii++)
        for(int jj = 0;jj<4;jj++)
        {
            E_matrix_0(ii,jj) = E[j][ii][jj];
            E_matrix_1(ii,jj) = E[j+1][ii][jj];

            A_matrix_0(ii,jj) = A[j][ii][jj];
            A_matrix_1(ii,jj) = A[j][ii][jj];
        }

    E_matrix_2.block(0,0,2,2) = E_matrix_1.block(0,0,2,2);
    E_matrix_2.block(0,2,2,2) = -E_matrix_0.block(0,2,2,2);
    E_matrix_2.block(2,0,2,2) = E_matrix_1.block(2,0,2,2);
    E_matrix_2.block(2,2,2,2) = -E_matrix_0.block(2,2,2,2);

    temp_matrix.block(0,0,2,2) = E_matrix_0.block(0,0,2,2)*A_matrix_0.block(0,0,2,2);
    temp_matrix.block(2,0,2,2) = E_matrix_0.block(2,0,2,2)*A_matrix_0.block(0,0,2,2);

    MatrixXcf temp1(4,4);
    temp1.block(0,0,4,4) = MatrixXcf::Zero(4,4);
    temp1.block(0,0,4,2) = E_matrix_2.inverse() * temp_matrix.block(0,0,4,2);
//    cout<<temp1<<endl;
    for(int ii = 0;ii<4;ii++)
        for(int jj = 0;jj<4;jj++)
        {
            this->RT[j][ii][jj] = temp1(ii,jj);
        }
}

void ScularFunction::calculate_A(float depth_z)
{
    float z = depth_z;
    for(int i = 1;i<media->layer - 1;i++)
    {
        this->A_d[i][0][0] = exp(-this->y[i]*(z - media->depth[i-1]));
        this->A_d[i][0][1] = 0.0;
        this->A_d[i][1][0] = 0.0;
        this->A_d[i][1][1] = exp(-this->v[i]*(z - media->depth[i-1]));

        this->A_u[i][0][0] = exp(-this->y[i]*(media->depth[i] - z));
        this->A_u[i][0][1] = 0.0;
        this->A_u[i][1][0] = 0.0;
        this->A_u[i][1][1] = exp(-this->v[i]*(media->depth[i] - z));
    }
    int i = media->layer - 1;
    this->A_d[i][0][0] = exp(-this->y[i]*(z - media->depth[i-1]));
    this->A_d[i][1][1] = exp(-this->v[i]*(z - media->depth[i-1]));
/*    for(int ii=0;ii<2;ii++){
        for(int jj=0;jj<2;jj++)
        {
            //fprintf(fp,"%6.5f  ",this->E[i][ii][jj]);
            cout<<this->A_d[i][ii][jj]<<endl;
        }
        //fprintf(fp,"\n");
    }
    cout<<endl;

    for(int ii=0;ii<2;ii++){
        for(int jj=0;jj<2;jj++)
        {
            //fprintf(fp,"%6.5f  ",this->E[i][ii][jj]);
            cout<<this->A_u[i][ii][jj]<<endl;
        }
        //fprintf(fp,"\n");
    }
    cout<<endl;*/
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
        this->v[i] = pow(complex<float>(pow(wavenumber,2) - pow(w/this->media->belta[i],2)),0.5);
        this->y[i] = pow(complex<float>(pow(wavenumber,2) - pow(w/this->media->alpha[i],2)),0.5);
    }
}

void ScularFunction::calculate_E(float frequency,float wavenumber)
{
    this->calculate_v_gamma(frequency,wavenumber);
    float w = 2*M_PI*frequency;
    complex<float>* X = this->malloc_array1_cf(this->media->layer);
    for(int i = 0;i<this->media->layer;i++)
    {
        X[i] = complex<float>(wavenumber*wavenumber) + pow(this->v[i],2);
    }
    FILE *fp = fopen("E.txt","w");
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

/*        for(int ii=0;ii<4;ii++){
            for(int jj=0;jj<4;jj++)
            {
                //fprintf(fp,"%6.5f  ",this->E[i][ii][jj]);
                cout<<this->E[i][ii][jj]<<endl;
            }
            //fprintf(fp,"\n");
        }
        int a = 1+2;*/
    }
    fclose(fp);

    this->release_array1_cf(X);
}

void ScularFunction::scular_output()
{
    for(int i = 0;i<this->media->N_k;i++)
    {
        fprintf(this->scular_out_fp,"%f %f \t",this->velocity[i],this->scular[i].real());
    }
    fprintf(this->scular_out_fp,"\n");
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
            for(int jj=0;jj<low;jj++)
                array3[i][ii][jj] = 0.0;
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
        for(int jj=0;jj<low;jj++)
            array2[i][jj] = 0.0;
    }
    return array2;
}

float* ScularFunction::malloc_array1(int low)
{
    float* array1 = (float*)malloc(sizeof(float)*low);
    for(int jj=0;jj<low;jj++)
        array1[jj] = 0.0;
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
            for(int jj=0;jj<low;jj++)
                array3[i][ii][jj] = complex<float>(0.0,0.0);
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
        for(int jj=0;jj<low;jj++)
            array2[i][jj] = complex<float>(0.0,0.0);
    }
    return array2;
}

complex<float>* ScularFunction::malloc_array1_cf(int low)
{
    complex<float>* array1 = (complex<float>*)malloc(sizeof(complex<float>)*low);
    for(int jj=0;jj<low;jj++)
        array1[jj] = complex<float>(0.0,0.0);
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
