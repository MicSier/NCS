#include "singlet.h"
#include "numericalmethods.h"
#include "cmath"
#include <iostream>
#include "tryplet.h"
#include "coupled.h"
#include <iomanip>
#include <random>

double coupled_fp(double x,double Tc,double mi,double gamma0,double tp)
{
    double ksip=ek(x,tp)-mi+g(x,gamma0),ksim=ek(x,tp)-mi-g(x,gamma0);
    double res,eps=0.000000000001;

    if(abs(ksim)<eps)
        res= 1.0/(2.0*Tc);
    else
        res= tanh(ksim/(2.0*Tc))/(2.0*ksim);

    if(abs(ksip)<eps)
        res+= 1.0/(2.0*Tc);
    else
        res+= tanh(ksip/(2.0*Tc))/(2.0*ksip);

    return res/2.0;
}

double coupled_fmdk(double x,double Tc,double mi,double gamma0,double tp)
{
    double ksip=ek(x,tp)-mi+g(x,gamma0),ksim=ek(x,tp)-mi-g(x,gamma0);
    double res,eps=0.000000000001;

    if(abs(ksim)<eps)
        res= 1.0/(2.0*Tc);
    else
        res= tanh(ksim/(2.0*Tc))/(2.0*ksim);

    if(abs(ksip)<eps)
        res-= 1.0/(2.0*Tc);
    else
        res-= tanh(ksip/(2.0*Tc))/(2.0*ksip);

    return -res/2.0*sqrt(d2(x));
}

double coupled_fpdk2 (double x,double Tc,double mi,double gamma0,double tp)
{
    double ksip=ek(x,tp)-mi+g(x,gamma0),ksim=ek(x,tp)-mi-g(x,gamma0);
    double res,eps=0.000000000001;

    if(abs(ksim)<eps)
        res= 1.0/(2.0*Tc);
    else
        res= tanh(ksim/(2.0*Tc))/(2.0*ksim);

    if(abs(ksip)<eps)
        res+= 1.0/(2.0*Tc);
    else
        res+= tanh(ksip/(2.0*Tc))/(2.0*ksip);

    return res/2.0*d2(x);
}

double coupled_gap(double Tc,double mi,double gamma0,int n,int k,double Vs,double Vt,double tp)
{
    const double pi=asin(1.0)*2.0;
    double a,b,c,temp;

    a=sc_integrate1D_gl_gpu(coupled_fp,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    b=sc_integrate1D_gl_gpu(coupled_fpdk2,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    c=sc_integrate1D_gl_gpu(coupled_fmdk,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    temp=(1.0-Vs*a)*(1.0-Vt*b)+Vs*Vt*c*c;

    return temp;

}

double coupled_deltas(double Tc,double mi,double gamma0,int n,int k,double Vs,double tp)
{
    const double pi=asin(1.0)*2.0;
    double a,b,deltaratiots;

    a=sc_integrate1D_gl_gpu(coupled_fp,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    b=sc_integrate1D_gl_gpu(coupled_fmdk,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    deltaratiots=(1.0-Vs*a)/(Vs*b);

    return deltaratiots;
}

double coupled_simplified_deltas(double Tc,double mi,double gamma0,int n,int k,double Vs,double Vt,double tp)
{
    const double pi=asin(1.0)*2.0;
    double a,b,c,deltaratiots;

    a=sc_integrate1D_gl_gpu(coupled_fpdk2,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    b=sc_integrate1D_gl_gpu(coupled_fmdk,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    c=sc_integrate1D_gl_gpu(coupled_fp,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    deltaratiots=Vt/Vs*a/b-Vt*c*a/b;

    return deltaratiots;
}

double coupled_deltat(double Tc,double mi,double gamma0,int n,int k,double Vt,double tp)
{
    const double pi=asin(1.0)*2.0;
    double a,b,deltaratiost;

    a=sc_integrate1D_gl_gpu(coupled_fpdk2,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    b=sc_integrate1D_gl_gpu(coupled_fmdk,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    deltaratiost=-(1.0-Vt*a)/(Vt*b);

    return deltaratiost;
}

double coupled_simplified_deltat(double Tc,double mi,double gamma0,int n,int k,double Vs,double Vt,double tp)
{
    const double pi=asin(1.0)*2.0;
    double a,b,c,deltaratiost;

    a=sc_integrate1D_gl_gpu(coupled_fp,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    b=sc_integrate1D_gl_gpu(coupled_fmdk,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    c=sc_integrate1D_gl_gpu(coupled_fpdk2,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
    deltaratiost=Vs/Vt*a/b-Vs*a*c/b;

    return deltaratiost;
}

double coupled_simplified_ratio(double Tcs,double mis,double Tct,double mit,double gamma0,int n,int k,double Vs,double Vt,double tp)
{
    const double pi=asin(1.0)*2.0;
    double fpd2s,fpd2t,fps,fpt,fmds,fmdt,ratio;

    fps=sc_integrate1D_gl_gpu(coupled_fp,0.0,pi,n,k,Tcs,mis,gamma0,tp)/(pi);
    fpt=sc_integrate1D_gl_gpu(coupled_fp,0.0,pi,n,k,Tct,mit,gamma0,tp)/(pi);
    fmds=sc_integrate1D_gl_gpu(coupled_fmdk,0.0,pi,n,k,Tcs,mis,gamma0,tp)/(pi);
    fmdt=sc_integrate1D_gl_gpu(coupled_fmdk,0.0,pi,n,k,Tct,mit,gamma0,tp)/(pi);
    fpd2s=sc_integrate1D_gl_gpu(coupled_fpdk2,0.0,pi,n,k,Tcs,mis,gamma0,tp)/(pi);
    fpd2t=sc_integrate1D_gl_gpu(coupled_fpdk2,0.0,pi,n,k,Tct,mit,gamma0,tp)/(pi);

  //  ratio=(Vt/Vs*fpd2s/fmds-Vt*(fps*fpd2s/fmds+fmds))*(Vs/Vt*fpt/fmdt-Vs*(fpt*fpd2t/fmdt+fmdt));
  //  ratio=fpd2s*fpt/(fmds*fmdt)-Vt*(fpd2s*fpt*fpd2t/(fmds*fmdt)+fpd2s*fmdt/fmds)-Vs*(fps*fpt*fpd2s/(fmds*fmdt)+fpt*fmds/fmdt)+Vs*Vt*(fps*fpd2s/fmds+fmds)*(fpt*fpd2t/fmdt+fmdt);
  //  cout<<fpd2s*fpt/(fmds*fmdt)<<" "<<-Vt*(fpd2s*fpt*fpd2t/(fmds*fmdt)+fpd2s*fmdt/fmds)<<" "<<-Vs*(fps*fpt*fpd2s/(fmds*fmdt)+fpt*fmds/fmdt)<<" "<<Vs*Vt*(fps*fpd2s/fmds+fmds)*(fpt*fpd2t/fmdt+fmdt)<<endl;
    ratio=(fpd2s*fpt-Vt*fpd2s*fpd2t*fpt-Vs*fps*fpt*fpd2s+Vs*Vt*fps*fpt*fpd2s*fpd2t)/(fmds*fmdt);
    return ratio;
}

double* coupled_get_res(double Vs,double Vt,double nl,double tp,double gamma0,double sxt,double syt,double sxc,double syc,int n,int k,bool isfromsinglet)
{
    double t,ch,chp,tol=0.00001;
    double xt,yt,x,y,at,bt,h;
    bool flag1=false;
    double* res_sin;
    double* res_try;
    double* res;
    int c=0;
    double ta[10];


    for(int i=0;i<10;i++)
    {
        ta[i]=0.0;
    }

    res=(double *)malloc(3*sizeof(double));

    res_sin=singlet_get_res(Vs,nl,tp,gamma0,sxt,syt,sxc,syc,n,k);
    res_try=tryplet_get_res(Vt,nl,tp,gamma0,sxt,syt,sxc,syc,n,k);

    xt=res_sin[0];
    yt=res_try[0];

    if(isfromsinglet)
    {
        at=xt;
        h=(yt-xt)/400.0;
    }
    else
    {
        at=yt;
        h=-(yt-xt)/400.0;
    }

     ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,at,gamma0,n,k,nl,tp);
     x=coupled_gap(at,ch,gamma0,n,k,Vs,Vt,tp);

      for(int i=1;i<=400 && flag1==false;i++)
      {
          bt=at+h;
          ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,bt,gamma0,n,k,nl,tp);
          y=coupled_gap(bt,ch,gamma0,n,k,Vs,Vt,tp);

          if(x*y<0)
          {
              if(isfromsinglet)
              {
                 xt=at;
                 yt=bt;
              }
              else
              {
                  xt=bt;
                  yt=at;
              }
              flag1=true;

          }

          x=y;
          at=bt;
      }

    t=coupled_solve1D_zbr(coupled_gap,xt,yt,tol,gamma0,n,k,Vs,Vt,sxc,syc,nl,tp);
    ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,n,k,nl,tp);
    chp=0.0;

    while(abs(ch-chp)>tol && t!=0.0 && c<10 && flag1)
    {
        chp=ch;
        t=coupled_solve1D_zbr(coupled_gap,xt,yt,tol,gamma0,n,k,Vs,Vt,sxc,syc,nl,tp);
        ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,n,k,nl,tp);
        for(int i=0;i<c;i++)
        {
            if(ta[i]==t)
             flag1=false;
        }
        ta[c]=t;
        c++;
    }

      res[0]=t;
      res[1]=ch;

      if(isfromsinglet)
        res[2]=coupled_deltas(t,ch,gamma0,n,k,Vs,tp);
      else
        res[2]=coupled_deltat(t,ch,gamma0,n,k,Vt,tp);


      std::cout<<"gamma0  "<<std::setprecision(12)<<gamma0<<"   ch   "<<std::setprecision(12)<<ch<<"  Tc  "<<std::setprecision(12)<<t<<"  dt/ds "<<std::setprecision(12)<<res[2]<<endl;

      return res;
}

double* eq_gamma0_fs(double gamma0,double Vs,double Vt,double nl,double tp,double sxc,double syc,int n,int k,double pts,double ptt,double pchs,double pcht,bool isfromsinglet)
{
    double* res;
    double dspdt,mi;
    double* tmp;

    tmp=(double *)malloc(4*sizeof(double));
    res=coupled_get_res_both(Vs,Vt,nl,tp,gamma0,sxc,syc,n,k,pts,ptt,pchs,pcht);

    tmp[0]=gamma0;
    if(isfromsinglet)
    {
        dspdt=1.0/res[3];
        mi=res[2];
    }
    else
    {
        dspdt=res[7];
        mi=res[6];
    }

    tmp[1]=-pow(mi,2.0)/4.0-pow(dspdt,2.0)+1.0+mi*gamma0/2.0*abs(dspdt)-pow(gamma0,2.0)/4.0*pow(dspdt,2.0);

    if(isfromsinglet)
    {
      tmp[2]=res[1];
      tmp[3]=res[2];
    }
    else
    {
      tmp[2]=res[5];
      tmp[3]=res[6];
    }

    return tmp;
}

double* eq_gamma0_bz(double gamma0,double Vs,double Vt,double nl,double tp,double sxc,double syc,int n,int k,double pts,double ptt,double pchs,double pcht,bool isfromsinglet)
{

    double* res;
    double dspdt;
    double* tmp;

    tmp=(double *)malloc(5*sizeof(double));
    res=coupled_get_res_both(Vs,Vt,nl,tp,gamma0,sxc,syc,n,k,pts,ptt,pchs,pcht);
    tmp[0]=gamma0;
    if(isfromsinglet)
    {
        dspdt=res[3];
    }
    else
    {
        dspdt=res[7];
    }

    tmp[1]=pow(dspdt,2.0)-1.0;

    if(isfromsinglet)
    {
      tmp[2]=res[1];
      tmp[3]=res[2];
      tmp[4]=res[3];
    }
    else
    {
      tmp[2]=res[5];
      tmp[3]=res[6];
      tmp[4]=res[7];
    }

cout<<"eq_gamma0_bz: "<<tmp[1]<<" t: "<<tmp[2]<<" mi: "<<tmp[3]<<" d: "<<tmp[4]<<endl;
    free(res);
    return tmp;
}

double* coupled_simplified_gap(double Tc,double mi,double gamma0,int n,int k,double Vs,double Vt,double tp)
{
  const double pi=asin(1.0)*2.0;
  double a,b,c,temp;
  double* res;

  res=(double *)malloc(6*sizeof(double));

  a=sc_integrate1D_gl_gpu(coupled_fp,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
  b=sc_integrate1D_gl_gpu(coupled_fpdk2,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);
  c=sc_integrate1D_gl_gpu(coupled_fmdk,0.0,pi,n,k,Tc,mi,gamma0,tp)/(pi);

  /*temp=1.0-Vs*a-Vt*b+Vs*Vt*c*c+Vs*Vt*a*b;
  cout<<-Vs*a<<" "<<-Vt*b<<" "<<Vs*Vt*c*c<<" "<<Vs*Vt*a*b<<" "<<temp<<endl;*/
  res[0]=-Vs*a;
  res[1]=-Vt*b;
  res[2]=Vs*Vt*a*b;
  res[3]=Vs*Vt*c*c;
  res[4]=1+res[0]+res[1]+res[2];
  res[5]=res[4]+res[3];

  return res;
}

double* coupled_get_res_both(double Vs,double Vt,double nl,double tp,double gamma0,double sxc,double syc,int n,int k,double pts,double ptt,double pchs,double pcht)
{
    double t,ch,chp,tol=0.0000001;
    double msxt,msyt,mtxt,mtyt,hs,ht,th,xs,ys,xt,yt,dt;
    bool flag1=true;
    double* res;
    int c=0,i=0;
    double ta[10];
    double sig;

    if(ptt>pts)
      sig=1.0;
    else
      sig=-1.0;

    res=(double *)malloc(8*sizeof(double));
    for(int i=0;i<10;i++)
    {
        ta[i]=0.0;
    }
    th=abs(ptt-pts)*0.1;
    t=0.0;
    dt=0.0;
    msxt=pts-sig*th;
    hs=th;
    double mh=(pcht-pchs)*0.001;

     ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,msxt,gamma0,n,k,nl,tp);
     xs=coupled_gap(msxt,ch,gamma0,n,k,Vs,Vt,tp);

    do
    {
        i++;
        msyt=msxt+sig*hs;
        ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,msyt,gamma0,n,k,nl,tp);
        ys=coupled_gap(msyt,ch,gamma0,n,k,Vs,Vt,tp);
        //cout<<"find s "<<msyt<<" "<<ys<<endl;
        if(xs*ys<0)
        {
            t=coupled_solve1D_zbr(coupled_gap,msxt,msyt,tol,gamma0,n,k,Vs,Vt,sxc,syc,nl,tp);
        }
        else
        {
            xs=ys;
            msxt=msyt;
        }

    } while (t==0.0 && i<100);

   if(t!=0.0)
   {
    ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,n,k,nl,tp);
   //cout<<"ch "<<ch<<endl;
    chp=0.0;

    while(abs(ch-chp)>tol && t!=0.0 && flag1)
    {
        chp=ch;
        t=coupled_solve1D_zbr(coupled_gap,msxt,msyt,tol,gamma0,n,k,Vs,Vt,sxc,syc,nl,tp);
        ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,n,k,nl,tp);

        for(int i=0;i<c;i++)
        {
            if(ta[i]==t)
             flag1=false;
        }
        ta[c]=t;
        c=(c+1)%10;
    }
    res[0]=gamma0;
    res[1]=t;
    res[2]=ch;
    res[3]=coupled_deltas(t,ch,gamma0,100*n,k,Vs,tp);
   }

else {
  cout<<i<<endl;
  res[0]=0.0;
  res[1]=0.0;
  res[2]=0.0;
  res[3]=0.0;
  return res;
}


      std::cout<<"singlet   :"<<std::setprecision(12)<<gamma0<<" "<<std::setprecision(12)<<t<<" "<<std::setprecision(12)<<ch<<" "<<std::setprecision(12)<<res[3]<<endl;

    for(int i=0;i<10;i++)
    {
        ta[i]=0.0;
    }
    c=0;
    t=0.0;
    i=0;
    ht=th*0.1;
    mtyt=ptt+sig*ht;
    //cout<<"sm0: "<<mtyt<<" "<<th<<" "<<ptt<<" "<<sig<<endl;
    ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,mtyt,gamma0,n,k,nl,tp);
    yt=coupled_gap(mtyt,ch,gamma0,n,k,Vs,Vt,tp);
    //cout<<"find t "<<mtyt<<" "<<yt<<endl;
    do
    {
        i++;
        mtxt=mtyt-sig*ht;

        ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,mtxt,gamma0,n,k,nl,tp);
        xt=coupled_gap(mtxt,ch,gamma0,n,k,Vs,Vt,tp);
        //cout<<"find t "<<mtxt<<" "<<xt<<endl;
        if(xt*yt<0.0)
        {
           //cout<<"sm0: "<<mtxt<<" "<<mtyt<<endl;
            t=coupled_solve1D_zbr(coupled_gap,mtxt,mtyt,tol,gamma0,n,k,Vs,Vt,sxc,syc,nl,tp);
        }
        else
        {
            yt=xt;
            mtyt=mtxt;
        }
        //cout<<"trip yt "<<yt<<" i "<<i<<endl;
    } while (t==0.0 && i<1000);

    if(t!=0.0)
    {
    ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,n,k,nl,tp);
    chp=0.0;
    flag1=true;

    while(abs(ch-chp)>tol && t!=0.0 && flag1)
    {
        chp=ch;
        t=coupled_solve1D_zbr(coupled_gap,mtxt,mtyt,tol,gamma0,n,k,Vs,Vt,sxc,syc,nl,tp);
        ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,n,k,nl,tp);
        //cout<<"mix_t t,ch "<<t<<" "<<ch<<endl;
        for(int i=0;i<c;i++)
        {
            if(ta[i]==t)
             flag1=false;
        }
        ta[c]=t;
        c=(c+1)%10;
    }

    res[4]=gamma0;
    res[5]=t;
    res[6]=ch;
    res[7]=coupled_deltat(t,ch,gamma0,100*n,k,Vt,tp);
    }
    else
    {
      cout<<i<<endl;
      res[3]=0.0;
      res[4]=0.0;
      res[5]=0.0;
    }


    std::cout<<"tryplet   :"<<std::setprecision(12)<<gamma0<<" "<<std::setprecision(12)<<t<<" "<<std::setprecision(12)<<ch<<" "<<std::setprecision(12)<<res[7]<<endl;
      return res;
}

/*double* coupled_get_res_both_simplified(double Vs,double Vt,double nl,double tp,double gamma0,double sxc,double syc,double xk[],double ak[],int n,int k,double pts,double ptt,double pchs,double pcht)
{
    double t,ch,chp,tol=0.0000001;
    double msxt,msyt,mtxt,mtyt,h,th,xs,ys,xt,yt,dt;
    bool flag1=true;
    double* res;
    int c=0,i=0;
    double ta[10];
    double sig;

    if(ptt>pts)
      sig=1.0;
    else
      sig=-1.0;

    res=(double *)malloc(6*sizeof(double));
    for(int i=0;i<10;i++)
    {
        ta[i]=0.0;
    }
    th=abs(ptt-pts)*0.05;
    t=0.0;
    dt=0.0;
    msxt=pts-sig*th;
    h=th;
    double mh=(pcht-pchs)*0.05;

     ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,msxt,gamma0,xk,ak,n,k,nl,tp);
     xs=coupled_simplified_gap(msxt,ch,gamma0,xk,ak,n,k,Vs,Vt,tp);

    do
    {
        i++;
        msyt=msxt+sig*h;
        ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,msyt,gamma0,xk,ak,n,k,nl,tp);
        ys=coupled_simplified_gap(msyt,ch,gamma0,xk,ak,n,k,Vs,Vt,tp);
        if(xs*ys<0)
        {
            t=coupled_solve1D_zbr(coupled_simplified_gap,msxt,msyt,tol,gamma0,xk,ak,n,k,Vs,Vt,sxc,syc,nl,tp);
            std::cout<<"t "<<t<<std::endl;
        }
        else
        {
            xs=ys;
            msxt=msyt;
        }

    } while (t==0.0 && i<100);

   if(t!=0.0)
   {
    ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,xk,ak,n,k,nl,tp);
   //cout<<"ch "<<ch<<endl;
    chp=0.0;

    while(abs(ch-chp)>tol && t!=0.0 && flag1)
    {
        chp=ch;
        t=coupled_solve1D_zbr(coupled_simplified_gap,msxt,msyt,tol,gamma0,xk,ak,n,k,Vs,Vt,sxc,syc,nl,tp);
        ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,xk,ak,n,k,nl,tp);

        for(int i=0;i<c;i++)
        {
            if(ta[i]==t)
             flag1=false;
        }
        ta[c]=t;
        c=(c+1)%10;
    }
   }


      res[0]=t;
      res[1]=ch;
      res[2]=coupled_deltas(t,ch,gamma0,xk,ak,n,k,Vs,tp);

      std::cout<<"singlet   :"<<std::setprecision(12)<<gamma0<<" "<<std::setprecision(12)<<t<<" "<<std::setprecision(12)<<ch<<" "<<std::setprecision(12)<<res[2]<<endl;

    for(int i=0;i<10;i++)
    {
        ta[i]=0.0;
    }
    c=0;
    t=0.0;
    i=0;
    mtyt=ptt+sig*th;
    ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,mtyt,gamma0,xk,ak,n,k,nl,tp);
    yt=coupled_simplified_gap(mtyt,ch,gamma0,xk,ak,n,k,Vs,Vt,tp);
    do
    {
        i++;
        mtxt=mtyt-sig*h;
        ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,mtxt,gamma0,xk,ak,n,k,nl,tp);
        xt=coupled_simplified_gap(mtxt,ch,gamma0,xk,ak,n,k,Vs,Vt,tp);
        if(xt*yt<0)
        {
          //  cout<<"sm0: "<<mtxt<<" "<<mtyt<<endl;
            t=coupled_solve1D_zbr(coupled_simplified_gap,mtxt,mtyt,tol,gamma0,xk,ak,n,k,Vs,Vt,sxc,syc,nl,tp);
        }
        else
        {
            yt=xt;
            mtyt=mtxt;
        }
    } while (t==0.0 && i<100);

    if(t!=0.0)
    {
    ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,xk,ak,n,k,nl,tp);
    chp=0.0;

    while(abs(ch-chp)>tol && t!=0.0 && flag1)
    {
        chp=ch;
        t=coupled_solve1D_zbr(coupled_simplified_gap,mtxt,mtyt,tol,gamma0,xk,ak,n,k,Vs,Vt,sxc,syc,nl,tp);
        ch=occ_solve1D_zbr(sc_occ,sxc,syc,tol,t,gamma0,xk,ak,n,k,nl,tp);
        //cout<<"mix_t t,ch "<<t<<" "<<ch<<endl;
        for(int i=0;i<c;i++)
        {
            if(ta[i]==t)
             flag1=false;
        }
        ta[c]=t;
        c=(c+1)%10;
    }

    }

    res[3]=t;
    res[4]=ch;
    res[5]=coupled_deltat(t,ch,gamma0,xk,ak,n,k,Vt,tp);

    std::cout<<"tryplet   :"<<std::setprecision(12)<<gamma0<<" "<<std::setprecision(12)<<t<<" "<<std::setprecision(12)<<ch<<" "<<std::setprecision(12)<<res[5]<<endl;
      return res;
}*/
