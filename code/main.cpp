#include <iostream>
#include "plotting.h"
#include "numericalmethods.h"
#include "singlet.h"
#include "tryplet.h"
#include "coupled.h"
#include <fstream>
#include <iomanip>
#include <sstream>

const double pi=3.141592653597932384;
string dn[] = {"coupled+s","coupled+t","univ+s","phase+","coupled+nodless+s","coupled+nodal+s","coupled+nodless+t","coupled+nodal+t","fs","nodes","univ+t","dos","simplification+gap","simplification+delta","simplification+ratio"};
double xk[]={0.9061798459386641,-0.9061798459386641,0.538469310105683,-0.5384693101056829,0.0};
double ak[]={0.2369268850561876,0.2369268850561876,0.47862867049936647,0.47862867049936586,0.5688888888888889};
double delta=0.0000000001;
int n=5000, k=5;

string tname (double Vs,double Vt,double nl,double tp)
{
    auto y = [] (double x)
    {
        ostringstream ss;
        ss<<x;

        return ss.str();
    };
    return "+V_s="+y(Vs)+"+V_t="+y(Vt)+"+nl="+y(nl)+"+t'="+y(tp);
}

string pname (double nl,double tp,double vmax)
{

    auto y = [] (double x)
    {
        ostringstream ss;
        ss<<x;

        return ss.str();
    };
    return  "+nl="+y(nl)+"+t'="+y(tp)+"+vmax="+y(vmax);
}

string uname (double ilr,double nl,double tp)
{
    auto y = [] (double x)
    {
        ostringstream ss;
        ss<<x;

        return ss.str();
    };
    return "+ilr="+y(ilr)+"+nl="+y(nl)+"+t'="+y(tp);
}

string dname (double tp,double g0)
{
    auto y = [] (double x)
    {
        ostringstream ss;
        ss<<x;

        return ss.str();
    };
    return "t'="+y(tp)+"+g_0'="+y(g0);
}

void check_univ(double nl,double ilr,double tp)
{
    double xt=0.0000000000001,yt=100.0,xc=-5.0,yc=5.0;
    double Vs,Vt,t,ch,chp,g,h;
    double* tabs;
    double* tabt;
    double* tab;
    double ptt,pts,pchs,pcht,vspvt=1.0;
    string names[3];
    int counter=0;


    for(int i=0;i<3;i++)
    {
        if(ilr<=1.0)
        {
            Vt=2.0-0.25*i;
            Vs=Vt*ilr;
        }
        else
        {
            Vs=2.0-0.25*i;
            Vt=Vs/ilr;
        }

        g=0.0;
        h=0.1;

        fstream os(path(dn[0]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream ot(path(dn[1]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream osn(path(dn[4]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream otn(path(dn[5]+tname(Vs,Vt,nl,tp),1),fstream::out);

        tabs=singlet_get_res(Vs,nl,tp,g,xt,yt,xc,yc,xk,ak,n,k);
        pts=tabs[0];
        pchs=tabs[1];
        tabt=tryplet_get_res(Vt,nl,tp,g,xt,yt,xc,yc,xk,ak,n,k);
        ptt=tabt[0];
        pcht=tabt[1];
        counter=0;
        cout<<"pts,ptt "<<pts<<" , "<<ptt<<endl;
        while(h>0.0000001)
        {
            g+=h;
            tab=coupled_get_res_both(Vs,Vt,nl,tp,g,xc,yc,xk,ak,n,k,pts,ptt,pchs,pcht);
            if(counter==10)
            {
                //h=h*10;
            }
            if(tab[0]<=0.0 || tab[3]<=0.0)
            {
                counter=0;
                g=g-h;
                cout<<"stare h "<<h<<" podaj nowe h"<<endl;
                cin>>h;
            }
            else
            {
                counter++;

                os<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<tab[0]<<" "<<std::setprecision(12)<<tab[1]<<" "<<std::setprecision(12)<<tab[2]<<endl;
                ot<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<tab[3]<<" "<<std::setprecision(12)<<tab[4]<<" "<<std::setprecision(12)<<tab[5]<<endl;

                pts=tab[0];
                ptt=tab[3];
                pchs=tab[1];
                pcht=tab[4];

            }



        }
        plot_tc(tname(Vs,Vt,nl,tp),dn[0]+tname(Vs,Vt,nl,tp),dn[1]+tname(Vs,Vt,nl,tp),dn[4]+tname(Vs,Vt,nl,tp),dn[5]+tname(Vs,Vt,nl,tp));
        names[i]=dn[2]+tname(Vs,Vt,nl,tp);
        coupled_universality_tabulate1D(dn[0]+tname(Vs,Vt,nl,tp),dn[4]+tname(Vs,Vt,nl,tp),dn[1]+tname(Vs,Vt,nl,tp),dn[5]+tname(Vs,Vt,nl,tp),names[i] );

    }

    plot_univ(dn[2]+"+dtpds"+uname(ilr,nl,tp),names[0],names[1],names[2],3);
    plot_univ(dn[6]+"+dspdt"+uname(ilr,nl,tp),names[0],names[1],names[2],2);
    //plot_univ(dn[2]+"+Tc"+uname(ilr,nl,tp),names[0],names[1],names[2],4);
    //plot_univ(dn[6]+"+Tc"+uname(ilr,nl,tp),names[0],names[1],names[2],5);


}

void check_simplification(double nl,double Vs,double Vt,double tp)
{
    double xt=0.0000000000001,yt=100.0,xc=0.000001,yc=10.0;
    double  t,ch,chp,g,h;
    double* tabs;
    double* tabt;
    double* tab1,*tab2;
    double* sgap,*tgap;
    double ptt,pts,pchs,pcht;
    double ptts,ptss,pchss,pchts;
    string names[3];
    int counter=0;

        g=0.0;
        h=0.1;

      /*  fstream os(path(dn[0]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream ot(path(dn[1]+tname(Vs,Vt,nl,tp),1),fstream::out);*/
        fstream sgs(path(dn[8]+tname(Vs,Vt,nl,tp)+"s",1),fstream::out);
        fstream sgt(path(dn[8]+tname(Vs,Vt,nl,tp)+"t",1),fstream::out);
        /*fstream sd(path(dn[9]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream sr(path(dn[10]+tname(Vs,Vt,nl,tp),1),fstream::out);*/
        fstream sc(path(dn[8]+tname(Vs,Vt,nl,tp)+"+cal+",1));

        tabs=singlet_get_res(Vs,nl,tp,g,xt,yt,xc,yc,xk,ak,n,k);
        pts=tabs[0];
        pchs=tabs[1];
        tabt=tryplet_get_res(Vt,nl,tp,g,xt,yt,xc,yc,xk,ak,n,k);
        ptt=tabt[0];
        pcht=tabt[1];
        counter=0;
        ptss=pts;
        ptts=ptt;
        pchss=pchs;
        pchts=pcht;

        while(h>0.0000001 && g<=0.2)
        {
            g+=h;
            tab1=coupled_get_res_both(Vs,Vt,nl,tp,g,xc,yc,xk,ak,n,k,pts,ptt,pchs,pcht);
            //tab1=coupled_get_res_both_simplified(Vs,Vt,nl,tp,g,xc,yc,xk,ak,n,k,ptss,ptts,pchss,pchts);

            if(counter==10)
            {
                //h=h*10;
            }
            if(tab1[0]<=0.0 || tab1[3]<=0.0)
            {
                counter=0;
                g=g-h;
                h=h/10.0;
            }
            else
            {
                counter++;

              /*  os<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<tab1[0]<<" "<<std::setprecision(12)<<tab1[1]<<" "<<std::setprecision(12)<<tab1[2]<<endl;
                ot<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<tab1[3]<<" "<<std::setprecision(12)<<tab1[4]<<" "<<std::setprecision(12)<<tab1[5]<<endl;
              */  sgap=coupled_simplified_gap(tab1[0],tab1[1],g,xk,ak,n,k,Vs,Vt,tp);
                tgap=coupled_simplified_gap(tab1[3],tab1[4],g,xk,ak,n,k,Vs,Vt,tp);
                sgs<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<sgap[0]<<" "<<std::setprecision(12)<<sgap[1]<<" "<<std::setprecision(12)<<sgap[2]<<" "<<std::setprecision(12)<<sgap[3]<<" "<<std::setprecision(12)<<sgap[4]<<" "<<std::setprecision(12)<<sgap[5]<<endl;
                sgt<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<tgap[0]<<" "<<std::setprecision(12)<<tgap[1]<<" "<<std::setprecision(12)<<tgap[2]<<" "<<std::setprecision(12)<<tgap[3]<<" "<<std::setprecision(12)<<tgap[4]<<" "<<std::setprecision(12)<<tgap[5]<<endl;
              /*  sd<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<tab1[2]<<" "<<std::setprecision(12)<<coupled_simplified_deltas(tab1[0],tab1[1],g,xk,ak,n,k,Vs,Vt,tp)<<" "<<std::setprecision(12)<<tab1[5]<<" "<<std::setprecision(12)<<coupled_simplified_deltat(tab1[3],tab1[4],g,xk,ak,n,k,Vs,Vt,tp)<<endl;
                sr<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<tab1[2]*tab1[5]<<" "<<std::setprecision(12)<<coupled_simplified_ratio(tab1[0], tab1[1], tab1[3], tab1[4],g,xk,ak,n,k, Vs, Vt, tp)<<endl;
              */  pts=tab1[0];
                ptt=tab1[3];
                pchs=tab1[1];
                pcht=tab1[4];
              /*  ptss=tab1[0];
                ptts=tab1[3];
                pchss=tab1[1];
                pchts=tab1[4];*/
                double fpd2s,fpd2t,fps,fpt,fmds,fmdt,ratio;

                fps=sc_integrate1D_gl_parallel(coupled_fp,0.0,pi,n,k,xk,ak,tab1[0],tab1[1],g,tp)/(pi);
                fpt=sc_integrate1D_gl_parallel(coupled_fp,0.0,pi,n,k,xk,ak,tab1[3],tab1[4],g,tp)/(pi);
                fmds=sc_integrate1D_gl_parallel(coupled_fmdk,0.0,pi,n,k,xk,ak,tab1[0],tab1[1],g,tp)/(pi);
                fmdt=sc_integrate1D_gl_parallel(coupled_fmdk,0.0,pi,n,k,xk,ak,tab1[3],tab1[4],g,tp)/(pi);
                fpd2s=sc_integrate1D_gl_parallel(coupled_fpdk2,0.0,pi,n,k,xk,ak,tab1[0],tab1[1],g,tp)/(pi);
                fpd2t=sc_integrate1D_gl_parallel(coupled_fpdk2,0.0,pi,n,k,xk,ak,tab1[3],tab1[4],g,tp)/(pi);
                //cout<<std::setprecision(12)<<g<<" "<<fps<<" "<<std::setprecision(12)<<fpt<<" "<<std::setprecision(12)<<fmds<<" "<<std::setprecision(12)<<fmdt<<" "<<std::setprecision(12)<<fpd2s<<" "<<std::setprecision(12)<<fpd2t<<endl;
                if(g==0.2)
                sc<<std::setprecision(12)<<Vs<<" "<<fps<<" "<<std::setprecision(12)<<fpt<<" "<<std::setprecision(12)<<fmds<<" "<<std::setprecision(12)<<fmdt<<" "<<std::setprecision(12)<<fpd2s<<" "<<std::setprecision(12)<<fpd2t<<endl;
            }



        }
        plot_tc(tname(Vs,Vt,nl,tp),dn[0]+tname(Vs,Vt,nl,tp),dn[1]+tname(Vs,Vt,nl,tp),dn[8]+tname(Vs,Vt,nl,tp),dn[9]+tname(Vs,Vt,nl,tp));
}

void univ_just_plot(double nl,double ilr,double tp)
{
  double xt=0.0000000000001,yt=100.0,xc=0.000001,yc=10.0;
  double Vs,Vt,t,ch,chp,g,h;
  double* tabs;
  double* tabt;
  double* tab;
  double ptt,pts,pchs,pcht,vspvt=1.0;
  string names[3];
  int counter=0;

  for(int i=0;i<3;i++)
  {
      if(ilr<=1.0)
      {
          Vt=2.0-0.5*i;
          Vs=Vt*ilr;
      }
      else
      {
          Vs=2.0-0.5*i;
          Vt=Vs/ilr;
      }
      names[i]=dn[2]+tname(Vs,Vt,nl,tp);
    }
  plot_univ(dn[2]+"+dtpds"+uname(ilr,nl,tp),names[0],names[1],names[2],3);
  plot_univ(dn[6]+"+dspdt"+uname(ilr,nl,tp),names[0],names[1],names[2],2);
}

void get_phase_diagram(double nl,double ilr,double tp,int np,double vmax,bool singlet)
{
    double xt=0.00000001,yt=10.0,xc=0.0,yc=5.0;
    double Vs,Vt,t,ch,chp,g,h,npbzs,npfss,npbzt,npfst,dspdt,dtpds;
    double g1,gh;
    double* tabs;
    double* tabt;
    double* tab;
    double eqfs=0.0,eqfsp=0.0,eqbs=0.0,eqbsp=0.0,eqft=0.0,eqftp=0.0,eqbt=0.0,eqbtp=0.0;
    double* gfs,*gbs,*gft,*gbt;
    double ptt,pts,pchs,pcht,vspvt=1.0;
    int count=0;
    string pha_names,pha_namet;
    bool fbs,fbt,ffs,fft;
    pha_names=dn[3]+pname(nl,tp,vmax)+"s";
    pha_namet=dn[3]+pname(nl,tp,vmax)+"t";

    fstream phas(path(pha_names,3),fstream::out);
    fstream phat(path(pha_namet,3),fstream::out);
    fstream phafs(path(pha_names+"fs",3),fstream::out);
    fstream phaft(path(pha_namet+"fs",3),fstream::out);



    for(int i=0;i<np;i++)
    {
      // check_univ(nl,ilr,tp);
        fbs=false;
        fbt=false;
        ffs=false;
        fft=false;
        if(singlet)
        {
          if(ilr<=1.0)
          {
              Vt=vmax;
              Vs=Vt*ilr;
          }
          else
          {
              Vs=vmax;
              Vt=Vs/ilr;
          }
        }
        else
        {
          if(ilr<=1.0)
          {
              Vs=vmax;
              Vt=Vs*ilr;
          }
          else
          {
              Vt=vmax;
              Vs=Vt/ilr;
          }
        }

        eqbsp=0.0;eqftp=0.0;eqbtp=0.0;eqbt=0.0;eqfsp=0.0;
        g=0.0;
        h=0.1;


        fstream os(path(dn[0]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream ot(path(dn[1]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream onodlesss(path(dn[4]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream onodals(path(dn[5]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream onodlesst(path(dn[6]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream onodalt(path(dn[7]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream ofs(path(dn[8]+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream onodess(path(dn[9]+"s"+tname(Vs,Vt,nl,tp),1),fstream::out);
        fstream onodest(path(dn[9]+"t"+tname(Vs,Vt,nl,tp),1),fstream::out);


        tabs=singlet_get_res(Vs,nl,tp,g,xt,yt,xc,yc,xk,ak,n,k);
        pts=tabs[0];
        pchs=tabs[1];
        tabt=tryplet_get_res(Vt,nl,tp,g,xt,yt,xc,yc,xk,ak,n,k);
        ptt=max(tabt[0],1e-8);
        pcht=tabt[1];
        free(tabs);
        free(tabt);
        cout<<"pts,ptt "<<pts<<" , "<<ptt<<endl;
        os<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<pts<<" "<<std::setprecision(12)<<pchs<<" "<<0.0<<endl;
        ot<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<ptt<<" "<<std::setprecision(12)<<pcht<<" "<<0.0<<endl;

        onodlesss<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<pts<<" "<<std::setprecision(12)<<pchs<<" "<<0.0<<endl;
        onodalt<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<ptt<<" "<<std::setprecision(12)<<pcht<<" "<<0.0<<endl;

        double fss1,fss2,fss3,fss4,fst1,fst2,fst3,fst4;
        fss1=fs(g,pchs,1.0);
        fss2=fs(g,pchs,-1.0);
        fss3=-fss1;
        fss4=-fss2;
        fst1=fs(g,pcht,1.0);
        fst2=fs(g,pcht,-1.0);
        fst3=-fss1;
        fst4=-fss2;
        ofs<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<fss1<<" "<<std::setprecision(12)<<fss2<<" "<<std::setprecision(12)<<fss3<<" "<<std::setprecision(12)<<fss4<<" "<<std::setprecision(12)<<fst1<<" "<<std::setprecision(12)<<fst2<<" "<<std::setprecision(12)<<fst3<<" "<<std::setprecision(12)<<fst4<<endl;

        double ns1,ns2,ns3,ns4,nt1,nt2,nt3,nt4;
        nt1=asin(0.0);
        nt2=-asin(0.0);
        nt3=pi-asin(0.0);
        nt4=-pi+asin(0.0);
        onodest<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<nt1<<" "<<std::setprecision(12)<<nt2<<" "<<std::setprecision(12)<<nt3<<" "<<std::setprecision(12)<<nt4<<endl;
        while(h>0.00000001)
        {

            g+=h;
            tab=coupled_get_res_both(Vs,Vt,nl,tp,g,xc,yc,xk,ak,n,k,pts,ptt,pchs,pcht);

            if(abs(tab[5]-tab[1])<=0.0000001 || tab[1]<=0.0 || tab[5]<=0.0)
            {
                g=g-h;
                /*cout<<"stare h "<<h<<" podaj nowe h"<<endl;
                cin>>h;*/
                h=h/10.0;
                count=0;
            }
            else
            {
              pts=tab[1];
              pchs=tab[2];
              dtpds=tab[3];

              ptt=tab[5];
              pcht=tab[6];
              dspdt=tab[7];

/*                if(count==20)
                {
                    h=h*10;
                    count=0;
                }*/

                count++;
                eqbs=pow(1.0/dtpds,2.0)-1.0;
                eqbt=pow(dspdt,2.0)-1.0;
                eqfs=-pow(pchs,2.0)/4.0-pow(1.0/dtpds,2.0)+1.0+pchs*g/2.0*abs(1.0/dtpds)-pow(g,2.0)/4.0*pow(1.0/dtpds,2.0);
                eqft=-pow(pcht,2.0)/4.0-pow(dspdt,2.0)+1.0+pcht*g/2.0*abs(dspdt)-pow(g,2.0)/4.0*pow(dspdt,2.0);

                //cout<<eqfs<<" "<<eqft<<endl;

              /*  if(eqfs*eqfsp<0.0 && !ffs)
                {
                    gfs=coupled_gamma0_solve1D_zbr(eq_gamma0_fs,g-h,g,delta,Vs,Vt,nl,tp,xc,yc,xk,ak,n,k,pts,ptt,pchs,pcht,true);
                    ffs=true;
                  //  onodals<<std::setprecision(12)<<gfs[0]<<" "<<std::setprecision(12)<<gfs[2]<<" "<<std::setprecision(12)<<<<" "<<std::setprecision(12)<<gfs[3]<<endl;
                  //  onodlesss<<std::setprecision(12)<<gfs[0]<<" "<<std::setprecision(12)<<gfs[2]<<" "<<std::setprecision(12)<<<<" "<<std::setprecision(12)<<gfs[3]<<endl;

                    //cout<<"if s"<<endl;
                }*/
                if(eqbs*eqbsp<0.0)
                {
                  /*  gbs=coupled_gamma0_solve1D_zbr(eq_gamma0_bz,g-h,g,eqbsp,eqbs,delta,Vs,Vt,nl,tp,xc,yc,xk,ak,n,k,pts,ptt,pchs,pcht,true);
                    fbs=true;
                    //phas<<ilr<<"  "<<setprecision(12)<<gbs[0]<<endl;
                    onodals<<std::setprecision(12)<<gbs[0]<<" "<<std::setprecision(12)<<gbs[2]<<" "<<std::setprecision(12)<<" "<<std::setprecision(12)<<gbs[3]<<" "<<std::setprecision(12)<<gbs[4]<<endl;
                    onodlesss<<std::setprecision(12)<<gbs[0]<<" "<<std::setprecision(12)<<gbs[2]<<" "<<std::setprecision(12)<<" "<<std::setprecision(12)<<gbs[3]<<" "<<std::setprecision(12)<<gbs[4]<<endl;

                    nt1=asin(1.0);
                    nt2=-asin(1.0);
                    nt3=pi-asin(1.0);
                    nt4=-pi+asin(1.0);
                    onodess<<std::setprecision(12)<<gbs[0]<<" "<<std::setprecision(12)<<nt1<<" "<<std::setprecision(12)<<nt2<<" "<<std::setprecision(12)<<nt3<<" "<<std::setprecision(12)<<nt4<<endl;
                    */
                    gh=h/100;
                    for(int i=0;i<100;i++)
                    {
                      g1=g-h+gh*i;
                      phas<<ilr<<"  "<<g1<<"  "<<setprecision(12)<<eq_gamma0_bz(g1,Vs,Vt,nl,tp,xc,yc,xk,ak,n,k,pts,ptt,pchs,pcht,true)[1]<<endl;
                    }

                }
              /*  if(eqft*eqftp<0.0 && !fft)
                {
                    gft=coupled_gamma0_solve1D_zbr(eq_gamma0_fs,g-h,g,delta,Vs,Vt,nl,tp,xc,yc,xk,ak,n,k,pts,ptt,pchs,pcht,false);
                    fft=true;
                    //  cout<<"if t"<<endl;
                }*/
                if(eqbt*eqbtp<0.0)
                {
                    /*gbt=coupled_gamma0_solve1D_zbr(eq_gamma0_bz,g-h,g,eqbtp,eqbt,delta,Vs,Vt,nl,tp,xc,yc,xk,ak,n,k,pts,ptt,pchs,pcht,false);
                    fbt=true;
                    //phat<<ilr<<"  "<<setprecision(12)<<gbt[0]<<endl;
                    onodalt<<std::setprecision(12)<<gbt[0]<<" "<<std::setprecision(12)<<gbt[2]<<" "<<std::setprecision(12)<<gbt[3]<<" "<<std::setprecision(12)<<gbt[4]<<endl;
                    onodlesst<<std::setprecision(12)<<gbt[0]<<" "<<std::setprecision(12)<<gbt[2]<<" "<<std::setprecision(12)<<gbt[3]<<" "<<std::setprecision(12)<<gbt[4]<<endl;

                    nt1=asin(1.0);
                    nt2=-asin(1.0);
                    nt3=pi-asin(1.0);
                    nt4=-pi+asin(1.0);
                    onodest<<std::setprecision(12)<<gbt[0]<<" "<<std::setprecision(12)<<nt1<<" "<<std::setprecision(12)<<nt2<<" "<<std::setprecision(12)<<nt3<<" "<<std::setprecision(12)<<nt4<<endl;
                    */
                    gh=h/100;
                    for(int i=0;i<100;i++)
                    {
                      g1=g-h+gh*i;
                      phat<<ilr<<"  "<<g1<<"  "<<setprecision(12)<<eq_gamma0_bz(g1,Vs,Vt,nl,tp,xc,yc,xk,ak,n,k,pts,ptt,pchs,pcht,false)[1]<<endl;
                    }

                }

                os<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<pts<<" "<<std::setprecision(12)<<pchs<<" "<<std::setprecision(12)<<dtpds<<endl;
                ot<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<ptt<<" "<<std::setprecision(12)<<pcht<<" "<<std::setprecision(12)<<dspdt<<endl;

                if(eqbs<0.0)
                {
                    onodals<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<pts<<" "<<std::setprecision(12)<<ptt<<" "<<std::setprecision(12)<<dtpds<<endl;
                }
                else
                {
                    onodlesss<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<pts<<" "<<std::setprecision(12)<<pchs<<" "<<std::setprecision(12)<<dtpds<<endl;
                }
                if(eqbt<0.0)
                {
                    onodalt<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<ptt<<" "<<std::setprecision(12)<<pcht<<" "<<std::setprecision(12)<<dspdt<<endl;

                }
                else
                {
                    onodlesst<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<ptt<<" "<<std::setprecision(12)<<pcht<<" "<<std::setprecision(12)<<dspdt<<endl;
                }

                fss1=fs(g,pchs,1.0);
                fss2=fs(g,pchs,-1.0);
                fss3=-fss1;
                fss4=-fss2;
                fst1=fs(g,pcht,1.0);
                fst2=fs(g,pcht,-1.0);
                fst3=-fst1;
                fst4=-fst2;
                ofs<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<fss1<<" "<<std::setprecision(12)<<fss2<<" "<<std::setprecision(12)<<fss3<<" "<<std::setprecision(12)<<fss4<<" "<<std::setprecision(12)<<fst1<<" "<<std::setprecision(12)<<fst2<<" "<<std::setprecision(12)<<fst3<<" "<<std::setprecision(12)<<fst4<<endl;

                if(abs(1.0/dtpds)<1)
                {
                  ns1=asin(abs(1.0/dtpds));
                  ns2=-asin(abs(1.0/dtpds));
                  ns3=pi-asin(abs(1.0/dtpds));
                  ns4=-pi+asin(abs(1.0/dtpds));
                  onodess<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<ns1<<" "<<std::setprecision(12)<<ns2<<" "<<std::setprecision(12)<<ns3<<" "<<std::setprecision(12)<<ns4<<endl;

                }
                if(abs(dspdt)<1)
                {
                  nt1=asin(abs(dspdt));
                  nt2=-asin(abs(dspdt));
                  nt3=pi-asin(abs(dspdt));
                  nt4=-pi+asin(abs(dspdt));
                  onodest<<std::setprecision(12)<<g<<" "<<std::setprecision(12)<<nt1<<" "<<std::setprecision(12)<<nt2<<" "<<std::setprecision(12)<<nt3<<" "<<std::setprecision(12)<<nt4<<endl;

                }



                eqfsp=eqfs;
                eqbsp=eqbs;
                eqftp=eqft;
                eqbtp=eqbt;
                free(tab);
            }


        }
            //phas<<ilr<<"  "<<g<<"  "<<setprecision(12)<<abs(tabs[0]-tabt[0])<<endl;
          /*  if(fbs)
            {
              phas<<ilr<<"  "<<setprecision(12)<<gbs[0]/g<<"  "<<setprecision(12)<<g<<"  "<<setprecision(12)<<gbs[0]<<"  "<<setprecision(12)<<gbs[1]<<"  "<<setprecision(12)<<abs(pts-ptt)<<endl;
              free(gbs);
            }
            if(fbt)
            {
              phat<<ilr<<"  "<<setprecision(12)<<gbt[0]/g<<"  "<<setprecision(12)<<g<<"  "<<setprecision(12)<<gbt[0]<<"  "<<setprecision(12)<<gbt[1]<<"  "<<setprecision(12)<<abs(pts-ptt)<<endl;
              free(gbt);
            }
            if(ffs)
              phafs<<ilr<<"  "<<setprecision(12)<<gfs[0]/g<<"  "<<setprecision(12)<<g<<"  "<<setprecision(12)<<gfs[0]<<"  "<<setprecision(12)<<gfs[1]<<endl;
            if(fft)
              phaft<<ilr<<"  "<<setprecision(12)<<gft[0]/g<<"  "<<setprecision(12)<<g<<"  "<<setprecision(12)<<gft[0]<<"  "<<setprecision(12)<<gft[1]<<endl;
*/
            phas<<ilr<<"  "<<setprecision(12)<<g<<"  "<<abs(pts-ptt)<<endl;
            phat<<ilr<<"  "<<setprecision(12)<<g<<"  "<<abs(pts-ptt)<<endl;

        //coupled_universality_tabulate1D(dn[0]+tname(Vs,Vt,nl,tp),dn[4]+tname(Vs,Vt,nl,tp),dn[1]+tname(Vs,Vt,nl,tp),dn[5]+tname(Vs,Vt,nl,tp),dn[2]+uname(ilr,nl,tp) );
        //plot_tc("Tc"+tname(Vs,Vt,nl,tp),dn[0]+tname(Vs,Vt,nl,tp),dn[4]+tname(Vs,Vt,nl,tp),dn[1]+tname(Vs,Vt,nl,tp),dn[5]+tname(Vs,Vt,nl,tp));

        ilr=ilr+0.01;
    }
  //  cout<<pha_names<<" "<<pha_namet<<endl;
      /*  plot_phase(pha_names,pha_names,true);
        plot_phase(pha_namet,pha_namet,false);
        /*plot_phase(pha_name+"fs",pha_name+"fs");
        plot_phase(pha_name+"ft",pha_name+"ft");*/
}

void phase_just_plot(double nl,double ilr,double tp,int np,double vmax)
{
  string pha_names,pha_namet;
  pha_names=dn[3]+pname(nl,tp,vmax)+"s";
  pha_namet=dn[3]+pname(nl,tp,vmax)+"t";
  plot_phase(pha_names,pha_names,true);
  plot_phase(pha_namet,pha_namet,false);
}

int main()
{
  double xt=0.000000001,yt=10.0,xc=-5.0,yc=5.0;
    double tp=0.0,nl=1.5,ilr=0.25,g0=0.01,ch,mi,dmi;
    string dos_name;
    double *res1,*res2;
    double ptt,pts,pchs,pcht;
    double* tabs;
    double* tabt;

    get_phase_diagram(1.5,0.25,tp,75,2,true);

    return 0;
}
