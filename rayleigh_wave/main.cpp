#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include"model.h"
#include"scularfunction.h"
#include "calculate_cofficient.h"
#include "dispersion_io.h"
#include "dispersioncurve.h"


int main()
{
    char* model_file = "/home/zgh/files/cuda/dispersion_curve/rayleigh_wave/input/model.txt";
    char* conf_file = "/home/zgh/files/cuda/dispersion_curve/rayleigh_wave/input/configure.txt";
    char* scular_file = "/home/zgh/files/cuda/dispersion_curve/rayleigh_wave/input/scular.txt";

    model* media = new model(model_file,conf_file);
    ScularFunction* scular = new ScularFunction(media,scular_file);
    dispersionCurve* dispersion = new dispersionCurve(scular);
    for(int i = 0 ;i<=media->N_f;i++)
    {
        cout<<"frequency is:"<<media->f[i]<<endl;
        dispersion->search_root(media->f[i]);        
        //scular->calculate_scular(media->f[i]);
        //scular->scular_output();
    }



    return 0;
}
