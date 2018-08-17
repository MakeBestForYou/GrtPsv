#ifndef MODEL_H
#define MODEL_H
#include<stdio.h>

class model
{
public:
    model(char* filename,char* conf_file);

public:

    char* model_file;
    char* conf_file;
    int layer;
    float* density;
    float* alpha;
    float* belta;
    float* depth;

    float* lambda;
    float* miu;
    float* xi;
    float* kxi;

    float* f;
    float pre_w;
    float pre_k;
    float* w;
    float* k;
    int N_k;
    int N_f;

    int N_mode;

    float maxVelocity;
    float minVelocity;

    float dVelocity;
    double deviation;

private:
    void model_read();
    void conf_read();

    void calculate_lambda();
    void calculate_miu();
    void calculate_xi_kxi();
};

#endif // MODEL_H
