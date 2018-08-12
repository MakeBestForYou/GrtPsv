#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"calculate_cofficient.h"

float* calculate_miu(float** model_par,int model_layer)
{
    float* miu = (float*)malloc(sizeof(float)*model_layer);
    for(int i = 0;i<model_layer;i++)
    {
        miu[i] = model_par[i][2]*model_par[i][2]*model_par[i][1]*pow(10,6);
    }
    return miu;
}

float* calculate_lambda(float** model_par,int model_layer)
{
    float* lambda = (float*)malloc(sizeof(float)*model_layer);
    for(int i = 0;i<model_layer;i++)
    {
        lambda[i] = model_par[i][3]*model_par[i][2]*model_par[i][1]*pow(10,6) - 2*model_par[i][2]*model_par[i][2]*model_par[i][1]*pow(10,6);
    }
    return lambda;
}

float* calculate_lambda(float** model_par,int model_layer);
