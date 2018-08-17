#include "dispersioncurve.h"

dispersionCurve::dispersionCurve(ScularFunction* scul)
{
    this->N_mode = scul->media->N_mode;
    this->scul = scul;
    curves= (float*)malloc(N_mode*sizeof(float));
    disper_fp = fopen("/home/zgh/files/cuda/dispersion_curve/rayleigh_wave/input/dispersion.txt","w+");
//    curves = (list<curve>*)malloc(sizeof(list<curve>));
}

void dispersionCurve::search_root(float fre)
{
    point present_point[3];
    for(int i=0;i<3;i++){
        present_point[i].x= -999;
        present_point[i].y = 1.0;
    }
    complex<float> tmp;
    float velocity = scul->media->minVelocity-0.5;

    bool flag = true;
    this->pre_mode = 0;

    fprintf(disper_fp,"%f   ",fre);
    while(flag)
    {
        float wavenumber = 2*M_PI*fre/velocity;

//        this->dVelocity = calculate_dv(fre,wavenumber);

        tmp = scul->calculate_scular(fre,wavenumber);
//        fprintf(fp,"%f  %f\n",velocity,tmp.real());
//        cout<<"phase velocity is "<<velocity<<" scul:"<<tmp.real()<<endl;
        refresh_data(present_point,tmp.real(),velocity);
        printf(" %f  %f %f     ",present_point[0].x,present_point[1].x,present_point[2].x);
        printf(" %f  %f %f\n\n\n",present_point[0].y,present_point[1].y,present_point[2].y);
        flag = try_search_root(&velocity,present_point,fre);
    }
    fprintf(disper_fp,"\n");
}

bool dispersionCurve::try_search_root(float *velocity,point* data,float frequency)
{
    //
    bool flag = true;

    float deviation = 0.00001;
    printf("%f\n",deviation);
    if(fabs(data[2].y)<deviation||fabs(data[1].y)<deviation||fabs(data[0].y)<deviation)
    {
        if(fabs(data[2].y)<deviation){
            addMode(frequency,data[2].x,this->pre_mode);
            printf("\n**********%f***************\n",data[2].y);
            data[2].y = 1;
        }else if(fabs(data[1].y)<deviation){
            addMode(frequency,data[1].x,this->pre_mode);
            printf("\n**********%f***************\n",data[1].y);
            data[1].y = -1;
        }else if(fabs(data[0].y)<deviation)
        {
            addMode(frequency,data[0].x,this->pre_mode);
            printf("\n**********%f***************\n",data[0].y);
            data[0].y = -1;
        }
        //addMode(frequency,*velocity,this->pre_mode);
        *velocity = data[2].x+scul->media->dVelocity;

        this->pre_mode++;
        printf("find one-------------------------------\n");
    }
    else
    {
        if(data[0].y>data[1].y&&(data[2].y>data[1].y))
        {
           *velocity = findSecondOrderMinPoint(data);
           //*velocity = *velocity+scul->media->dVelocity;
        }
        else
        {
            *velocity = *velocity+scul->media->dVelocity;
            //*velocity = *velocity+this->dVelocity;
        }
    }

    if(*velocity>scul->media->maxVelocity+1.0||this->pre_mode>(this->N_mode-1))
        flag = false;
    return flag;
}

float dispersionCurve::findSecondOrderMinPoint(point* data)
{
    float min_x;
    MatrixXd A(3,3);
    MatrixXd x(3,1);
    MatrixXd y(3,1);

    A<<double(data[0].x)*double(data[0].x),double(data[0].x),1.0,
       double(data[1].x)*double(data[1].x),double(data[1].x),1.0,
       double(data[2].x)*double(data[2].x),double(data[2].x),1.0;
    y<<double(data[0].y),double(data[1].y),double(data[2].y);
    x=A.inverse()*y;
    min_x = float(-x(1)/(2*x(0)));

//    cout<<x<<endl;
//    cout<<A<<endl;
//    cout<<A.inverse()<<endl;

//    printf("%f  %f %f :%f\n",data[0].x,data[1].x,data[2].x,min_x);
//    printf("%f  %f %f :%f\n",data[0].y,data[1].y,data[2].y,min_x);
//    cout<<min_x<<endl;
//    exit(1);

    if(min_x<data[2].x)
        return min_x;
    else
    {
        printf("find second order min x error!");
        exit(1);
        return min_x;
    }
}

void dispersionCurve::addMode(float f, float c, int Mode)
{
//    curve curve_tmp;
//    curve_tmp.c = c;
//    curve_tmp.f = f;

    curves[Mode] = c;
    fprintf(disper_fp,"%f   ",c);
    //curves[Mode].push_back(curve_tmp);
}

void dispersionCurve::refresh_data(point* data_pos,float sculf_value,float velocity)
{
    if(velocity<data_pos[2].x)
    {
        if(velocity<data_pos[1].x)
        {
            data_pos[2] = data_pos[1];
            data_pos[1].x = velocity;
            data_pos[1].y = sculf_value;
        }else
        {
            data_pos[2].x = velocity;
            data_pos[2].y = sculf_value;
        }
    }else
    {
        data_pos[0] = data_pos[1];
        data_pos[1] = data_pos[2];
        data_pos[2].x = velocity;
        data_pos[2].y = sculf_value;
    }
}

float dispersionCurve::calculate_dv(float frequency,float wavenumber)
{
    float k = wavenumber;
    float f = frequency;
    float w = 2*M_PI*f;

    complex<float>* v = (complex<float>*)malloc(sizeof(complex<float>)*scul->media->layer);
    complex<float>* y = (complex<float>*)malloc(sizeof(complex<float>)*scul->media->layer);

    for(int ii=0;ii<scul->media->layer;ii++)
    {
        v[ii] = pow(complex<float>(k*k-(w/scul->media->alpha[ii])*(w/scul->media->alpha[ii])),0.5);
        y[ii] = pow(complex<float>(k*k-(w/scul->media->belta[ii])*(w/scul->media->belta[ii])),0.5);
//        cout<<v[ii]<<endl;
//        cout<<y[ii]<<endl;
    }

    float N=0;
    for(int ii=1;ii<scul->media->layer;ii++)
    {
        N = N+2*(v[ii].imag()+y[ii].imag());
    }
    N = N/scul->media->belta[1];

    float dc = (scul->media->belta[2]-scul->media->minVelocity)/N;
//    cout<<N<<endl;
    return dc;

}
