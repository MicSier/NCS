#include "singlet.h"
#include "numericalmethods.h"
#include "cmath"
#include <iostream>
#include "tryplet.h"

double gd(double x,double gamma0)
{
    return gamma0*(sin(x)*sin(x));
  //  return gamma0*abs(sin(x));
}

double d2(double x)
{
    return pow(sin(x),2.0);
}

double tryplet_uifunt(double x,double Tc,double mi,double gamma0,double tp)
{
    double ksi=ek(x,tp)-mi,eps=0.000000001,vg=g(x,gamma0);
    double uifunt1=0.0,uifunt2=0.0,pare=0.0,uifunt;
  //  cout<<endl;
      if((abs(ksi)<eps)&&(abs(vg)<eps))
      {
                   uifunt1=1.0/(2.0*Tc);
                   //cout<<"uifunt 1"<<endl;
      }
      else
      {
          if(abs(ksi)<eps)
          {
              uifunt1=tanh(vg/(2.0*Tc))/vg;
              pare=d2(x)-pow(gd(x,gamma0)/vg,2.0);

              if(abs(vg/Tc)>100)
              {
                  uifunt2=-1.0/vg;
                  // cout<<"uifunt 2"<<endl;
              }
              else if(abs(vg/Tc)<eps)
                    {
                        uifunt2=0.0;
                //         cout<<"uifunt 3"<<endl;
                    }
                    else
                    {
                      uifunt2=1.0/(Tc+Tc*cosh(vg/Tc))-sinh(vg/Tc)/(vg+vg*cosh(vg/Tc));
              //         cout<<"uifunt 4"<<endl;
                    }

          }
          else
          {
              if(abs(vg)<eps)
              {
                  uifunt1=tanh((ksi)/(2.0*Tc))/(ksi);
            //      cout<<"uifunt 5"<<endl;
              }
              else if((abs(ksi+vg)<eps)||(abs(ksi-vg)<eps))
                    {
                        uifunt1=(Tc*tanh(vg/Tc)+vg)/(4.0*Tc*vg);
                        uifunt2=0.250*(tanh(vg/Tc)/vg-1.0/Tc);
                        pare=d2(x)-pow(gd(x,gamma0)/vg,2.0);
          //              cout<<"uifunt 6"<<endl;
                    }
                    else
                    {
                        uifunt1=tanh((ksi+vg)/(2.0*Tc))/(2.0*(ksi+vg))+ tanh((ksi-vg)/(2.0*Tc))/(2.0*(ksi-vg));
                        uifunt2=vg/ksi*(tanh((ksi+vg)/(2.0*Tc))/(2.0*(ksi+vg))-tanh((ksi-vg)/(2.0*Tc))/(2.0*(ksi-vg)));
                        pare=d2(x)-pow(gd(x,gamma0)/vg,2.0);
        //                cout<<"uifunt 7"<<endl;
                    }

          }

      }
      uifunt=uifunt1*d2(x)+uifunt2*pare;
      ///cout<<"uifunt "<<uifunt<<" uifunt1 "<<uifunt1<<" uifunt2 "<<uifunt2<<" pare "<<pare<<endl;
      return uifunt;
}

double tryplet_gap(double Tc,double mi,double gamma0,int n,int k,double Vt,double tp)
{
    const double pi=asin(1.0)*2.0;
    double temp;

    temp=sc_integrate1D_gl_gpu(tryplet_uifunt,0.0,pi,n,k,Tc,mi,gamma0,tp);
    temp=1.0-Vt/(2.0*pi)*(temp);
    return temp;

}

double* tryplet_get_res(double Vt,double nl,double tp,double gamma0,double xt,double yt,double xc,double yc,int n,int k)
{
    double t,ch,chp,tol=0.000000001;
    double* res;
  res=sc_solve1D_zbr(tryplet_gap,xt,yt,tol,gamma0,n,k,Vt,xc,yc,nl,tp);

    /*cout<<"t tryplet ch "<<chp<<" t "<<t<<endl;
    t=sc_solve1D_zbr(tryplet_gap,xt,yt,tol,0.0,gamma0,n,k,Vt);
    cout<<"ch tryplet ch "<<chp<<" t "<<t<<endl;
    ch=sc_solve1D_zbr(sc_occ,xc,yc,tol,t,gamma0,n,k,nl);

    chp=0.0;

    while(abs(ch-chp)>tol)
    {
        chp=ch;
        cout<<"t tryplet ch "<<chp<<" t "<<t<<endl;
        t=sc_solve1D_zbr(tryplet_gap,xt,yt,tol,chp,gamma0,n,k,Vt);
        cout<<"ch tryplet ch "<<chp<<" t "<<t<<endl;
        ch=sc_solve1D_zbr(sc_occ,xc,yc,tol,t,gamma0,n,k,nl);
    }

        res[0]=t;
        res[1]=ch;*/
      return res;
}
