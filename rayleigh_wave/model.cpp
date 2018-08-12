#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include "model.h"
#include "dispersion_io.h"

model::model(char* model_file,char* conf_file)
{
    this->model_file = model_file;
    this->model_read();
    this->calculate_miu();
    this->calculate_lambda();

    this->conf_file = conf_file;
    this->conf_read();
}

void model::model_read()
{
    FILE* fp;
    if(!(fp = fopen(this->model_file,"r")))
    {
        printf("model file open error");
        exit(1);
    }

    nextline(fp,4);
    fscanf(fp,"%d",&(this->layer));
    nextline(fp,2);

    this->alpha = (float*)(malloc(sizeof(float)*this->layer));
    this->belta = (float*)(malloc(sizeof(float)*this->layer));
    this->depth = (float*)(malloc(sizeof(float)*this->layer));
    this->density = (float*)(malloc(sizeof(float)*this->layer));

    int temp;
    for(int i = 0;i<this->layer;i++)
    {
        fscanf(fp,"%d",&temp);
        fscanf(fp,"%f",&(this->depth[i]));
        fscanf(fp,"%f",&(this->density[i]));
        fscanf(fp,"%f",&(this->belta[i]));
        fscanf(fp,"%f",&(this->alpha[i]));
        nextline(fp,1);
    }

    fclose(fp);
}

void model::conf_read()
{
    FILE* fp;
    if(!(fp = fopen(this->conf_file,"r")))
    {
        printf("configure file open error");
        exit(1);
    }
    float start_f;
    float end_f;
    float df;
    float start_k;
    float end_k;
    float dk;

    nextline(fp,1);
    fscanf(fp,"%f",&start_f);
    fscanf(fp,"%f",&end_f);
    fscanf(fp,"%f",&df);
    nextline(fp,2);
    fscanf(fp,"%f",&start_k);
    fscanf(fp,"%f",&end_k);
    fscanf(fp,"%f",&dk);

    this->N_f = int(fabs(start_f-end_f)/df+0.5);
    this->f = (float*)(malloc(sizeof(float)*this->N_f));
    this->w = (float*)(malloc(sizeof(float)*this->N_f));

    this->N_k = int(fabs(start_k-end_k)/dk+0.5);
    this->k = (float*)(malloc(sizeof(float)*this->N_k));

    for(int i = 0;i<=this->N_f;i++)
    {
        this->f[i] = start_f+i*df;
        this->w[i] = 2*M_PI*this->f[i];
    }
    for(int i = 0;i<this->N_k;i++)
    {
        this->k[i] = start_k+i*dk;
    }

    fclose(fp);
}

void model::calculate_miu()
{
    this->miu = (float*)malloc(sizeof(float)*this->layer);
    for(int i = 0;i<this->layer;i++)
    {
        this->miu[i] = (this->alpha[i])*(this->alpha[i])*(this->density[i]);
    }
}

void model::calculate_lambda()
{
    this->lambda = (float*)malloc(sizeof(float)*this->layer);
    for(int i = 0;i<this->layer;i++)
    {
        this->lambda[i] = this->belta[i]*this->belta[i]*this->density[i] - 2*this->miu[i];
    }
}

void model::calculate_xi_kxi()
{
    this->xi = (float*)malloc(sizeof(float)*this->layer);
    this->kxi = (float*)malloc(sizeof(float)*this->layer);
    for(int i = 0;i<this->layer;i++)
    {
        this->xi[i] = 1/(this->lambda[i]+2*this->miu[i]);
        this->kxi[i] = 4*this->miu[i]*(this->lambda[i]+this->miu[i])/(this->lambda[i]+2*this->miu[i]);
    }
}
