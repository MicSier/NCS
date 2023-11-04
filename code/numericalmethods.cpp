#include "numericalmethods.h"
#include "plotting.h"
#include <chrono>
#include <numeric>
#include <vector>
#include <execution>
#include "singlet.h"

double occ_solve1D_zbr(double (*f)(double,double,double,double*,double*,int,int,double,double),double a, double b,double tol,double mi,double gamma0,double xk[],double ak[],int n,int k,double Vs,double tp)
{
    int itmax=100;
    double d,r,s,e,p,q,xm,tol1,c,fa,fb,fc,eps=3.0e-8;

    fa=f(a,mi,gamma0,xk,ak,n,k,Vs,tp);
    fb=f(b,mi,gamma0,xk,ak,n,k,Vs,tp);

    if(fa*fb>0.0)
    {
        //cout<<"occ_zbr err:Takie same znaki!  fa  "<<fa<<" fb  "<<fb<<" tc,g0,a,b  "<<mi<<"  "<<gamma0<<"  "<<a<<"  "<<b<<endl;
        return 0.0;
    }
    c=b;
    fc=fb;
    for(int i=1;i<itmax;i++)
    {
        if(fb*fc>0.0)
        {
            c=a;
            fc=fa;
            d=b-a;
            e=d;
        }
        if(abs(fc)<abs(fb))
        {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*eps*abs(b)+0.5*tol;
        xm=0.5*(c-b);
        if((abs(xm)<tol1)||(fb==0)) return b;
        if((abs(e)>tol1)&&(abs(fa)>abs(fb)))
        {
            s=fb/fa;
            if(a==c)
            {
                p=2.0*xm*s;
                q=1.0-s;
            }
            else
            {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if(p>0.0) q=-q;
            p=abs(p);
            if(2.0*p<min(3.0*xm*q-abs(tol1*q),abs(e*q)))
            {
                e=d;
                d=p/q;
            }
            else
            {
                d=xm;
                e=d;
            }
        }
        else
        {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if(abs(d)>tol1) b+=d;
        else
        {
            if(xm>0.0) b=b+abs(tol1);
            else b=b-abs(tol1);
        }
        fb=f(b,mi,gamma0,xk,ak,n,k,Vs,tp);


    }
    cout<<"zbr exeding max iteractions!"<<endl;
    return b;
}

double* sc_solve1D_zbr(double (*f)(double,double,double,double*,double*,int,int,double,double),double a, double b,double tol,double gamma0,double xk[],double ak[],int n,int k,double Vs,double sxc,double syc,double nl,double tp)
{
    int itmax=100;
    double d,r,s,e,p,q,xm,tol1,c,fa,fb,fc,eps=3.0e-8,mi;
    double* res;

    res=(double *)malloc(2*sizeof(double));

    mi=occ_solve1D_zbr(sc_occ,sxc,syc,tol,a,gamma0,xk,ak,n,k,nl,tp);
    fa=f(a,mi,gamma0,xk,ak,n,k,Vs,tp);
     mi=occ_solve1D_zbr(sc_occ,sxc,syc,tol,b,gamma0,xk,ak,n,k,nl,tp);
    fb=f(b,mi,gamma0,xk,ak,n,k,Vs,tp);
    res[0]=b;
    res[1]=mi;

    if(fa*fb>0.0)
    {
        cout<<"sc_zbr err:Takie same znaki!  fa  "<<fa<<" fb  "<<fb<<" mi,g0,a,b  "<<mi<<"  "<<gamma0<<"  "<<a<<"  "<<b<<endl;
        res[0]=0.0;
        res[1]=0.0;
        return res;
    }
    c=b;
    fc=fb;
    for(int i=1;i<itmax;i++)
    {
        if(fb*fc>0.0)
        {
            c=a;
            fc=fa;
            d=b-a;
            e=d;
        }
        if(abs(fc)<abs(fb))
        {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*eps*abs(b)+0.5*tol;
        xm=0.5*(c-b);
        if((abs(xm)<tol1)||(fb==0)) return res;
        if((abs(e)>tol1)&&(abs(fa)>abs(fb)))
        {
            s=fb/fa;
            if(a==c)
            {
                p=2.0*xm*s;
                q=1.0-s;
            }
            else
            {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if(p>0.0) q=-q;
            p=abs(p);
            if(2.0*p<min(3.0*xm*q-abs(tol1*q),abs(e*q)))
            {
                e=d;
                d=p/q;
            }
            else
            {
                d=xm;
                e=d;
            }
        }
        else
        {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if(abs(d)>tol1) b+=d;
        else
        {
            if(xm>0.0) b=b+abs(tol1);
            else b=b-abs(tol1);
        }
         mi=occ_solve1D_zbr(sc_occ,sxc,syc,tol,b,gamma0,xk,ak,n,k,nl,tp);
        fb=f(b,mi,gamma0,xk,ak,n,k,Vs,tp);
        res[0]=b;
        res[1]=mi;

    }
    cout<<"zbr exeding max iteractions!"<<endl;
    return res;
}

double coupled_solve1D_zbr(double (*f)(double,double,double,double*,double*,int,int,double,double,double),double a, double b,double tol,double gamma0,double xk[],double ak[],int n,int k,double Vs,double Vt,double sxc,double syc,double nl,double tp)
{
    int itmax=100;
    double mi,d,r,s,e,p,q,xm,tol1,c,fa,fb,fc,eps=3.0e-8;

    mi=occ_solve1D_zbr(sc_occ,sxc,syc,tol,a,gamma0,xk,ak,n,k,nl,tp);
    fa=f(a,mi,gamma0,xk,ak,n,k,Vs,Vt,tp);
    mi=occ_solve1D_zbr(sc_occ,sxc,syc,tol,b,gamma0,xk,ak,n,k,nl,tp);
    fb=f(b,mi,gamma0,xk,ak,n,k,Vs,Vt,tp);


    if(fa*fb>0.0)
    {
        if(min(fa,fb)<tol*tol) return min(fa,fb);
        else
        {
        cout<<"coupled_zbr err:Takie same znaki!  fa  "<<fa<<"  fb  "<<fb<<endl;
        return 0.0;
        }

    }
    c=b;
    fc=fb;
    for(int i=1;i<itmax;i++)
    {
        if(fb*fc>0.0)
        {
            c=a;
            fc=fa;
            d=b-a;
            e=d;
        }
        if(abs(fc)<abs(fb))
        {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*eps*abs(b)+0.5*tol;
        xm=0.5*(c-b);
        if((abs(xm)<tol1)||(fb==0)) return b;
        if((abs(e)>tol1)&&(abs(fa)>abs(fb)))
        {
            s=fb/fa;
            if(a==c)
            {
                p=2.0*xm*s;
                q=1.0-s;
            }
            else
            {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if(p>0.0) q=-q;
            p=abs(p);
            if(2.0*p<min(3.0*xm*q-abs(tol1*q),abs(e*q)))
            {
                e=d;
                d=p/q;
            }
            else
            {
                d=xm;
                e=d;
            }
        }
        else
        {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if(abs(d)>tol1) b+=d;
        else
        {
            if(xm>0.0) b+=abs(tol1);
            else b+=-abs(tol1);
        }
        mi=occ_solve1D_zbr(sc_occ,sxc,syc,tol,b,gamma0,xk,ak,n,k,nl,tp);
        fb=f(b,mi,gamma0,xk,ak,n,k,Vs,Vt,tp);

    }
    cout<<"zbr exeding max iteractions!"<<endl;
    return b;
}

double* coupled_gamma0_solve1D_zbr(double* (*f)(double,double,double,double,double,double,double,double*,double*,int,int,double,double,double,double,bool),double a, double b,double fgmh,double fg,double tol,double Vs,double Vt,double nl,double tp,double sxc,double syc,double xk[],double ak[],int n,int k,double pts,double ptt,double pchs,double pcht,bool isfromsinglet)
{
    int itmax=100;
    double d,r,s,e,p,q,xm,tol1,c,fa,fb,fc,eps=3.0e-8;
    double* taba,*tabb,*tabc;
    //taba=f(a,Vs,Vt,nl,tp,sxc,syc,xk,ak,n,k,pts,ptt,pchs,pcht,isfromsinglet);
    fa=fgmh;
  //  tabb=f(b,Vs,Vt,nl,tp,sxc,syc,xk,ak,n,k,pts,ptt,pchs,pcht,isfromsinglet);
    fb=fg;

    if(fa*fb>0.0)
    {
        cout<<"gammacrit_zbr err:Takie same znaki! a "<<a<<" b "<<b<<" fa  "<<fa<<"  fb  "<<fb<<endl;
        tabb[0]=0.0;
        tabb[1]=0.0;
        return tabb;

    }
    c=b;
    fc=fb;
    //tabc=tabb;
    for(int i=1;i<itmax;i++)
    {
        if(fb*fc>0.0)
        {
            c=a;
            fc=fa;
            tabc=taba;
            d=b-a;
            e=d;
        }
        if(abs(fc)<abs(fb))
        {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
            taba=tabb;
            tabb=tabc;
            tabc=taba;
        }
        tol1=2.0*eps*abs(b)+0.5*tol;
        xm=0.5*(c-b);
        if((abs(xm)<tol1)||(fb==0)) return tabb;
        if((abs(e)>tol1)&&(abs(fa)>abs(fb)))
        {
            s=fb/fa;
            if(a=c)
            {
                p=2.0*xm*s;
                q=1.0-s;
            }
            else
            {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if(p>0.0) q=-q;
            p=abs(p);
            if(2.0*p<min(3.0*xm*q-abs(tol1*q),abs(e*q)))
            {
                e=d;
                d=p/q;
            }
            else
            {
                d=xm;
                e=d;
            }
        }
        else
        {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if(abs(d)>tol1) b+=d;
        else
        {
            if(xm>0.0) b+=abs(tol1);
            else b+=-abs(tol1);
        }
        tabb=f(b,Vs,Vt,nl,tp,sxc,syc,xk,ak,n,k,pts,ptt,pchs,pcht,isfromsinglet);
        fb=tabb[1];
    }
    cout<<"zbr exeding max iteractions!"<<endl;
    return tabb;
}


double sc_solve1D_new(double (*f)(double,double,double,double*,double*,int,int,double),double a,double tol,double mi,double gamma0,double xk[],double ak[],int n,int k,double Vs)
{
    double x,h,dev,vy;

    auto y = [f,mi,gamma0,Vs,xk,ak,n,k] (double t)
    {
        return f(t,mi,gamma0,xk,ak,n,k,Vs);
    };

    auto deriv = [y] (double t,double dt)
    {
        return (y(t+dt)-y(t-dt))/(2.0*dt);
    };

    x=a;
    dev=deriv(x,tol*10.0);
    vy=y(x);
    h=vy/dev;

    while(abs(h)>=tol)
    {
        h=vy/dev;
        x=x-h;

        vy=y(x);
        dev=deriv(x,tol);
        std::cout<<" h "<<h<<"   x   "<<x<<"  y  "<<vy<<std::endl;
    }

    return x;

}

double coupled_solve1D_new(double (*f)(double,double,double,double*,double*,int,int,double,double),double a,double tol,double mi,double gamma0,double xk[],double ak[],int n,int k,double Vs,double Vt)
{
    double x,h,dev,vy;
    int i=0;

    auto y = [f,mi,gamma0,Vs,Vt,xk,ak,n,k] (double t)
    {
        return f(t,mi,gamma0,xk,ak,n,k,Vs,Vt);
    };

    auto deriv = [y] (double t,double dt)
    {
        return (y(t+dt)-y(t-dt))/(2.0*dt);
    };

    x=a;
    dev=deriv(x,tol);
    vy=y(x);
    h=vy/dev;

    while(abs(h)>=tol && i<50)
    {
        h=vy/dev*0.2;
        x=x-h;

        vy=y(x);
        dev=deriv(x,tol);
        //std::cout<<" h "<<h<<"   x   "<<x<<"  y  "<<vy<<" i "<<i<<std::endl;
        i++;
    }

    if(i>=49)
    {
        return 0.0;
    }
    else
    {
       return x;
    }

}

double sc_integrate1D_gl_parallel(double (*f)(double,double,double,double,double),double a, double b, int n, int k,double* xk,double* ak,double Tc,double mi,double gamma0,double tp)
{
  std::vector<double> v(n*k);
   double h=(b-a)/(n);

  auto y = [a,h,xk,ak,f,Tc,mi,gamma0,tp] (int i)
{
  int j=i%5,m=i/5;
  double z=((2.0*a+h*(2.0*m+1.0))-h*xk[j])*0.5;
    return ak[j]*f(z,Tc,mi,gamma0,tp);
};

  std::iota(v.begin(), v.end(), 0);

  /*std::transform(v.begin(),v.end(),v.begin(),y);
  const double result = std::reduce(std::execution::par,v.begin(),v.end())*h*0.5;*/
  const double result = std::transform_reduce(std::execution::par,v.begin(), v.end(),0.0,[](auto a, auto b) {return a + b;},y)*h*0.5;

  return result;

}
