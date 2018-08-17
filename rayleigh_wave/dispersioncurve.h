#ifndef DISPERSIONCURVE_H
#define DISPERSIONCURVE_H
#include "scularfunction.h"
#include <list>

using namespace std;


struct curve
{
    float f;
    float c;
};

struct point
{
    float x;
    float y;
};

class dispersionCurve
{
public:
    dispersionCurve(ScularFunction* scul);
    void search_root(float fre);
public:
    int N_mode;
    int pre_mode;
    float* curves;
    //list<curve>* curves;
    model* media;
    ScularFunction* scul;
    float dVelocity;
    FILE* disper_fp;

private:
    void addMode(float f,float c,int Mode);

    void refresh_data(point* data_pos,float sculf_value,float velocity);
    bool try_search_root(float* wavenumber,point* tmp,float frequency);
    float findSecondOrderMinPoint(point* data);
    float calculate_dv(float frequency,float wavenumber);
};

#endif // DISPERSIONCURVE_H
