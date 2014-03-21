#include <iostream>
#include <cstdio>
#include "vec_ops.h"
#include <math.h>
using namespace std;


const double G = 6.67e-11;
const double step = 1;
const double time = 10;
const int bodies = 2;
double m[bodies];
VECTOR a[bodies];
VECTOR v[bodies];
VECTOR r[bodies];
VECTOR rcap[bodies][bodies];
double rdif[bodies][bodies];

FILE *fa,*fv,*fr;


void input_vector(VECTOR v[],int bodies)
{
    for(int i=0;i<bodies;i++)
        cin>>v[i].vec[0]>>v[i].vec[1]>>v[i].vec[2];
}


void input_scalar(double s[],int bodies)
{
    for(int i=0;i<bodies;i++)
        cin>>s[i];
}


void calculate_rdif_rcap()
{
    for(int i=0;i<bodies;i++)
    {
        for(int j=0;j<bodies-1;j++)
        {
            rdif[i][j] = sqrt( pow(r[j].vec[0]-r[i].vec[0],2) + pow(r[j].vec[1]-r[i].vec[1],2) + pow(r[j].vec[1]-r[i].vec[1],2) );
            rcap[i][j].vec[0] = (r[j].vec[0]-r[i].vec[0]) / rdif[i][j];
            rcap[i][j].vec[1] = (r[j].vec[1]-r[i].vec[1]) / rdif[i][j];
            rcap[i][j].vec[2] = (r[j].vec[2]-r[i].vec[2]) / rdif[i][j];
        }
    }

}


VECTOR new_vector(VECTOR x,VECTOR y)
{
    VECTOR t;
    for(int i=0;i<3;i++)
    {
        t.vec[i]=x.vec[i]+step*y.vec[i];
    }
    return t;
}

VECTOR equation_solver(VECTOR ac,int j)
{
    double Gxm_j=G*m[j];
    for(int i=0;i<bodies;i++)
    {
        if(i==j) continue;
        for(int k=0;k<3;k++)
        {
            ac.vec[k]= ac.vec[k]+ ((Gxm_j/rdif[i][j])*rcap[i][j].vec[k]);
        }
    }
    return ac;
}

VECTOR rk4_acceleration(VECTOR ac,int j)
{
    VECTOR k1,k2,k3,k4,t;
    int i;
    k1=equation_solver(ac,j);
    for(i=0;i<3;i++) t.vec[i]=ac.vec[i]+step*k1.vec[i]/2;
    k2=equation_solver(t,j);
    for(i=0;i<3;i++) t.vec[i]=ac.vec[i]+step*k2.vec[i]/2;
    k3=equation_solver(t,j);
    for(i=0;i<3;i++) t.vec[i]=ac.vec[i]+step*k3.vec[i];
    k4=equation_solver(t,j);
    for(i=0;i<3;i++) ac.vec[i] = ac.vec[i] + step/6*(k1.vec[i]+2*k2.vec[i]+2*k3.vec[i]+k4.vec[i]);
    cout<<ac.vec[0]<<ac.vec[1]<<ac.vec[2]<<endl;
    return ac;
}



void simulate_bodies()
{
    for(double i=0;i<time;i=i+step)
    {
        for(int j=0;j<bodies;j++)
        {
        a[j]=rk4_acceleration(a[j],j);
        v[j]=new_vector(v[j],a[j]);
        r[j]=new_vector(r[j],v[j]);
        fprintf(fa,"%f\t%f\t%f\t%f\n",i,a[j].vec[0],a[j].vec[1],a[j].vec[2]);
        fprintf(fv,"%f\t%f\t%f\t%f\n",i,v[j].vec[0],v[j].vec[1],v[j].vec[2]);
        fprintf(fr,"%f\t%f\t%f\t%f\n",i,r[j].vec[0],r[j].vec[1],r[j].vec[2]);
        }
        calculate_rdif_rcap();
    }
}

int main()
{
    cout<<"input mass"<<endl;
    input_scalar(m,bodies);
    cout<<"input acceleration"<<endl;
    input_vector(a,bodies);
    cout<<"input velocity"<<endl;
    input_vector(v,bodies);
    cout<<"input position"<<endl;
    input_vector(r,bodies);
    fa=fopen("a_vector.txt","w");
    fv=fopen("v_vector.txt","w");
    fr=fopen("r_vector.txt","w");
    fprintf(fa,"Time\ti\tj\tk\n");
    fprintf(fv,"Time\ti\tj\tk\n");
    fprintf(fr,"Time\ti\tj\tk\n");
    calculate_rdif_rcap();
    simulate_bodies();
}
