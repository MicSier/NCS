#include "singlet.h"
#include "numericalmethods.h"
#include "cmath"
#include <iostream>

const double pi=3.141592653597932384;

/*double evth(double x)
{
    if(abs(x)>100)
        return x/abs(x);
    else
    {
        double z=exp(x);
        return (z-1.0/z)/(z+1.0/z);
    }
}*/

double g(double x,double gamma0)
{
    return gamma0*(abs(sin(x)));
    //return gamma0*1;
}

double ek(double x,double tp)
{
    return -2.0*(cos(x))+2.0*tp*cos(2*x);
}

__device__ double sc_uifunc(double x,double Tc,double mi,double gamma0,double tp)
{
    double ksip=ek(x,tp)-mi+g(x,gamma0),ksim=ek(x,tp)-mi-g(x,gamma0);

    return tanh(ksip/(2.0*Tc))+tanh(ksim/(2.0*Tc));
}

__device__ double singlet_uifunt(double x,double Tc,double mi,double gamma0,double tp)
{
    double ksip=ek(x,tp)-mi+g(x,gamma0),ksim=ek(x,tp)-mi-g(x,gamma0);
    double res;

    if(abs(ksim)<0.0000000000001)
        res= 1.0/(2.0*Tc);
    else
        res= tanh(ksim/(2.0*Tc))/(2.0*ksim);

    if(abs(ksip)<0.0000000000001)
        res+= 1.0/(2.0*Tc);
    else
        res+= tanh(ksip/(2.0*Tc))/(2.0*ksip);

    return res;
}

 double singlet_gap(double Tc,double mi,double gamma0,int n,int k,double Vs,double tp)
{
    const double pi=asin(1.0)*2.0;
    double temp;

    temp=sc_integrate1D_gl_gpu(singlet_uifunt,0.0,pi,n,k,Tc,mi,gamma0,tp);
    temp=1.0-Vs/(2.0*pi)*(temp);
    return temp;

}

double sc_occ(double mi,double Tc,double gamma0,int n,int k,double nl,double tp)
{
    double temp;
    const double pi=3.141592653597932384;

    temp=sc_integrate1D_gl_gpu(sc_uifunc,0.0,pi,n,k,Tc,mi,gamma0,tp);

    temp=1.0-nl-temp/(2.0*pi);

    return temp;
}

double* singlet_get_res(double Vs,double nl,double tp,double gamma0,double xt,double yt,double xc,double yc,int n,int k)
{
    double t,ch,chp,tol=0.000000001;
    double* res;

    std::cout<<"begin singlet_get_res"<<std::endl;
    res=sc_solve1D_zbr(singlet_gap,xt,yt,tol,gamma0,n,k,Vs,xc,yc,nl,tp);
    std::cout<<"end singlet_get_res"<<std::endl;
    std::cout<<"t singlet ch "<<res[0]<<" t "<<res[1]<<std::endl;
   /* chp=0.0;
    cout<<"t singlet ch "<<chp<<" t "<<t<<endl;
    t=sc_solve1D_zbr(singlet_gap,xt,yt,tol,chp,gamma0,n,k,Vs);
    cout<<"ch singlet ch "<<chp<<" t "<<t<<endl;
    ch=sc_solve1D_zbr(sc_occ,xc,yc,tol,t,gamma0,n,k,nl);

    while(abs(ch-chp)>tol)
    {
        chp=ch;
        cout<<"t singlet ch "<<chp<<" t "<<t<<endl;
        t=sc_solve1D_zbr(singlet_gap,xt,yt,tol,chp,gamma0,n,k,Vs);
        cout<<"ch singlet ch "<<chp<<" t "<<t<<endl;
        ch=sc_solve1D_zbr(sc_occ,xc,yc,tol,t,gamma0,n,k,nl);

    }

      res[0]= t;
      res[1]=ch;*/

      return res;
}

double* sc_dos(double o,double g0,double tp)
{
  double* res;
  int npi=100000,nomega=100000;
  double ddos=1.e-02;
	double step=2.0*pi/npi;
	double densitytot=2.*npi;
  double omega,gamma,omp,omm;

  res=(double *)malloc(3*sizeof(double));


	double domega=2.*ddos;

    double densityp=0.;
    double densitym=0.;

    double x=-pi;
    for(int j=0;j<npi;j++)
    {
      omega=ek(x,tp);
      gamma=g(x,g0);
      omp=omega+gamma;
      omm=omega-gamma;
      if (abs(omp-o) < ddos) densityp=densityp+1.;
      if (abs(omm-o) < ddos) densitym=densitym+1.;
      x=x+step;
    }
    densitym=densitym/densitytot;
    densityp=densityp/densitytot;
    double densitypm=densityp+densitym;
    densitym=densitym/domega;
    densityp=densityp/domega;
    densitypm=densitypm/domega;

    res[0]=densitym;
    res[1]=densityp;
    res[2]=densitypm;

    return res;
}

double fs(double gamma0,double mi,double lambda)
{
  double a=1+0.25*gamma0*gamma0-0.25*mi*mi;
  if(a>0)
  {
    a=-0.5*mi+lambda*0.5*gamma0*sqrt(a);
  }
  else
  {
    a=-0.5*mi;
  }
  a=a/(1+0.25*gamma0*gamma0);
  return a;
}
