#include "plotting.h"
#include "singlet.h"
#include <iomanip>
#include <cmath>

const double pi=3.141592653597932384;

string path(string name,int n)
{
  string so="rashba";
    switch(n)
    {
        case 0:
        return "\"./data/"+so+"/tc/"+name+".dat\"";

        case 1:
        return "./data/"+so+"/tc/"+name+".dat";

        case 2:
        return "\"./data/"+so+"/phase/"+name+".dat\"";

        case 3:
        return "./data/"+so+"/phase/"+name+".dat";

        case 4:
        return "\"./data/"+so+"/uni/"+name+".dat\"";

        case 5:
        return "./data/"+so+"/uni/"+name+".dat";

        case 6:
        return "\"./data/"+so+"/dos/"+name+".dat\"";

        case 7:
        return "./data/"+so+"/dos/"+name+".dat";

        case 8:
        return "\"./plots/"+so+"/tc/"+name+".png\"";

        case 9:
        return "\"./plots/"+so+"/phase/"+name+".eps\"";

        case 10:
        return "\"./plots/"+so+"/uni/"+name+".png\"";

        case 11:
        return "\"./plots/"+so+"/dos/"+name+".eps\"";

        default:
        return "";
    }
}

void tabulate1D(string fname,double (*f)(double),double a, double b, int N)
{
    fstream outfile(path(fname,2),fstream::out);
    double x;
    double h=(b-a)/(N-1);

    if(!outfile.good())
    {
        cout<<"dobrze otwarty plik!"<<endl;
    }

    for(int i=0;i<N;i++)
    {
        x=a+i*h;
        outfile<<x<<" "<<f(x)<<endl;
    }
    outfile.close();
}

void plot_tc(string nameout,string namein1,string namein2,string namein3,string namein4)
{
    gnuplot p;

    p("set term png");
    p("set output "+path(nameout,8));
    p("set grid");
    p("set xlabel \"gamma_0\"");
    //p("set xrange [$3:$4]");
    p("set ylabel \"T_c\"");
    p("set title \""+nameout+"\"");
    //p("plot "+path(namein1,0)+" u 1:2 w l lw 2 title \""+namein1+"\", "+path(namein2,0)+" u 1:2 w l lw 2 title \""+namein2+"\" ,"+path(namein3,0)+" u 1:2 w l lw 2 title \""+namein3+"\" , "+path(namein4,0)+" u 1:2 w l lw 2 title \""+namein4+"\"");
    p("plot "+path(namein1,0)+" u 1:2 w l lw 2 title \""+namein1+"\", "+path(namein2,0)+" u 1:2 w l lw 2 title \""+namein2+"\" ,"+path(namein3,0)+" u 1:2 w p title \""+namein3+"\" , "+path(namein4,0)+" u 1:2 w p title \""+namein4+"\"");

}

void plot_phase(string nameout,string namein,bool isfromsinglet)
{
    gnuplot p;

    p("set term eps");
    p("set output "+path(nameout,9));
    p("set xlabel \"V_s/V_t\"");
    //p("set xrange [$3:$4]");
    p("set ylabel \"g_c/g_m\"");
    p("set title \""+nameout+"\"");
    p("set style fill transparent solid 0.2");
    p("set grid");
    p("set linetype 1");
    p("set key outside");

    if(isfromsinglet)
      p("plot "+path(namein,2)+" w l lw 3 title \"critical line\", "+path(namein,2)+" using ($1):($2) with filledcurves x1 title \"not nodal\", "+path(namein,2)+" using ($1):($2):(1) with filledcurves lt 7 title \"nodal\" ");
    else
      p("plot "+path(namein,2)+" w l lw 3 title \"critical line\", "+path(namein,2)+" using ($1):($2) with filledcurves x1 lt 7 title \"nodal\", "+path(namein,2)+" using ($1):($2):(1) with filledcurves title \"not nodal\" ");

}

void plot_univ(string nameout,string namein1,string namein2,string namein3,int n)
{
    gnuplot p;
    string numb=to_string(n);
    //string yl[]={"\"d_t/d_s\"","\"d_s/d_t\"","\"T_c/T_0\"","\"T_c/T_0\""};

    p("set term png");
    p("set output "+path(nameout,10));
    p("set grid");
    //p("set ylabel a");//+yl[n]
    //p("set xrange [$3:$4]");
    p("set xlabel \"g_0/g_m\"");
    p("set title \""+nameout+"\"");
    p("plot "+path(namein1,4)+" u 1:"+numb+" w p title \""+namein1+"\", "+path(namein2,4)+" u 1:"+numb+" w p title \""+namein2+"\","+path(namein3,4)+" u 1:"+numb+" w p title \""+namein3+"\"");
    /*p("plot "+path(namein2,0)+" u 1:"+numb+" w p");
    p("plot "+path(namein3,0)+" u 1:"+numb+" w p "); //title \""+namein3+"\"*/
}

void plot_dos(string nameout,string namein,double mi,double V)
{
    gnuplot p;
    ostringstream ssmi,ssvl,ssvp;
    ssmi<<mi;
    ssvl<<mi-V;
    ssvp<<mi+V;
    p("set term eps");
    p("set output "+path(nameout,11));
    p("set grid");
    p("set xlabel \"omega\"");
    //p("set xrange [$3:$4]");
    p("set ylabel \"dos\"");
    p("set style fill transparent solid 0.3");
    p("set title \""+nameout+"\"");
    p("filter(x,max) = (x < max) ? x : 1/0");
    p("plot "+path(namein,6)+" u 1:4 w l title \"density of states\", \'\' using (filter($1,"+ssmi.str()+")):($4) with filledcurves x1 lt 4 title \"occupied states\", \'-\' w p lt 7 title \"mi\", \'\' w p lt 6 title \"range\"");
    p(ssmi.str()+" 0.0");
    p("e");
    p(ssvl.str()+" 0.0");
    p(ssvp.str()+" 0.0");
    p("e");

}

void sc_tabulate1D(string fname,double* (*f)(double,double,double,double,double,double,double,double,double*,double*,int,int),double a, double b, int N,double Vs,double tp,double nl,double xt,double yt,double xc,double yc,double xk[],double ak[],int n,int k)
{
    fstream outfile(path(fname,2),fstream::out);
    double g0=a;
    double h=(b-a)/(N-1);
    double* tmp;

    if(!outfile.good())
    {
        cout<<"nie otwarty plik!"<<endl;
    }

    while(h>0.0001)
    {
        g0=g0+h;
        tmp=f(Vs,nl,tp,g0,xt,yt,xc,yc,xk,ak,n,k);
        if(tmp[0]==0.0)
        {
            g0=g0-h;
            h=h/2.0;
        }
        else
        {
            outfile<<setprecision(10)<<g0<<" "<<setprecision(10)<<tmp[0]<<" "<<setprecision(10)<<tmp[1]<<endl;
        }

    }
  }

void coupled_tabulate1D(string fname,double* (*f)(double,double,double,double,double,double,double,double,double*,double*,int,int,bool),double a, double b, int N,double Vs,double Vt,double nl,double xt,double yt,double xc,double yc,double xk[],double ak[],int n,int k,bool isfromsinglet)
{
    fstream outfile(path(fname,1),fstream::out);
    double g0=a;
    double h=(b-a)/(N-1);
    double* tmp;

    if(!outfile.good())
    {
        cout<<"nie otwarty plik!"<<endl;
    }

    while(h>0.000000001)
    {
        g0=g0+h;
        tmp=f(Vs,Vt,nl,g0,xt,yt,xc,yc,xk,ak,n,k,isfromsinglet);
        if(tmp[0]==0.0)
        {
            g0=g0-h;
            h=h/2.01;
        }
        else
        {
            outfile<<std::setprecision(12)<<g0<<" "<<std::setprecision(12)<<tmp[0]<<" "<<std::setprecision(12)<<tmp[1]<<" "<<std::setprecision(12)<<tmp[2]<<endl;
        }

    }

    outfile.close();
}

void coupled_gamma0_fs_tabulate1D(string fname,double* (*f)(double,double,double,double,double,double,double,double*,double*,int,int,bool,double,double),double a, double b, int N,double Vs,double Vt,double nl,double xt,double yt,double xc,double yc,double xk[],double ak[],int n,int k,bool isfromsinglet)
{
    fstream outfile(path(fname,2),fstream::out);
    double g0=a;
    double h=(b-a)/(N-1);
    double* tmp;

    if(!outfile.good())
    {
        cout<<"nie otwarty plik!"<<endl;
    }

    while(h>0.0000001)
    {
        g0=g0+h;
        tmp=f(Vs,Vt,nl,xt,yt,xc,yc,xk,ak,n,k,isfromsinglet,a,b);

        if(tmp[0]==0.0)
        {
            g0=g0-h;
            h=h/2.0;
        }
        else
        {
            outfile<<std::setprecision(12)<<g0<<" "<<std::setprecision(12)<<tmp[0]<<" "<<std::setprecision(12)<<tmp[1]<<" "<<std::setprecision(12)<<tmp[2]<<endl;
        }

    }

    outfile.close();
}

void coupled_universality_tabulate1D(string fname_singlet1,string fname_singlet2,string fname_tryplet1,string fname_tryplet2, string fname)
{
    fstream infile_s1(path(fname_singlet1,1));
    fstream infile_s2(path(fname_singlet2,1));
    fstream infile_t1(path(fname_tryplet1,1));
    fstream infile_t2(path(fname_tryplet2,1));

    fstream outfile(path(fname,5),fstream::out);


    string word;
    double g0_s[400],dtpds[400],g0_t[400],dspdt[400];
    double t_s[400],t_t[400];
    double mi_s[400],mi_t[400];
    int i=0,j=0,minindi,minindj,m,k;
    double minv,temp;

    while(infile_s1.good() && i<400)
    {
        infile_s1>>word;
        g0_s[i]=stod(word);
        infile_s1>>word;
        t_s[i]=stod(word);
        infile_s1>>word;
        mi_s[i]=stod(word);
        infile_s1>>word;
        dtpds[i]=stod(word);
        if(infile_s1.good()) i++;
    }

    while(infile_s2.good() && i<400)
    {
        infile_s2>>word;
        g0_s[i]=stod(word);
        infile_s2>>word;
        t_s[i]=stod(word);
        infile_s2>>word;
        mi_s[i]=stod(word);
        infile_s2>>word;
        dtpds[i]=stod(word);
        if(infile_s2.good()) i++;
    }

    while(infile_t1.good() && j<400)
    {
        infile_t1>>word;
        g0_t[j]=stod(word);
        infile_t1>>word;
        t_t[j]=stod(word);
        infile_t1>>word;
        mi_t[j]=stod(word);
        infile_t1>>word;
        dspdt[j]=stod(word);
       if(infile_t1.good()) j++;
    }

    while(infile_t2.good() && j<400)
    {
        infile_t2>>word;
        g0_t[j]=stod(word);
        infile_t2>>word;
        t_t[j]=stod(word);
        infile_t2>>word;
        mi_t[j]=stod(word);
        infile_t2>>word;
        dspdt[j]=stod(word);
       if(infile_t2.good()) j++;
    }

    minindi=0;
    minindj=0;
    minv=abs(dspdt[minindj]*dtpds[minindi]-1.0);

  /*  for(k=0;k<=i;k++)
    {*/
        for(m=1;m<j;m++)
        {
            temp=abs(dspdt[m]*dtpds[m]-1.0);
            if(temp<minv)
            {
                minindi=m;
                minindj=m;
                minv=temp;
            }

        }
  //  }

    double gammacrit=g0_t[i-1];
    cout<<i<<" "<<j<<" "<<minindi<<" "<<minindj<<endl;
    cout<<"gc  "<<gammacrit<<"  "<<g0_s[minindi]<<"  "<<g0_t[minindj]<<" "<<endl;
    for(k=0;k<i;k++)
    {
        outfile<<setprecision(10)<<g0_s[k]/gammacrit<<" "<<setprecision(10)<<pow(dspdt[k],2.0)+4*pow(mi_t[k],2.0)<<" "<<setprecision(10)<<pow(1.0/dtpds[k],2.0)+4*pow(mi_t[k],2.0)<<" "<<endl;//setprecision(10)<<t_s[k]/t_s[0]<<" "<<setprecision(10)<<t_t[k]/t_t[0]<<" "<<setprecision(10)<<mi_s[k]/mi_s[0]<<" "<<setprecision(10)<<mi_t[k]/mi_t[0]<<endl;
        //outfile<<setprecision(10)<<g0_s[k]/gammacrit<<" "<<setprecision(10)<<dspdt[k]<<" "<<setprecision(10)<<dtpds[k]<<" "<<endl;//setprecision(10)<<t_s[k]/t_s[0]<<" "<<setprecision(10)<<t_t[k]/t_t[0]<<" "<<setprecision(10)<<mi_s[k]/mi_s[0]<<" "<<setprecision(10)<<mi_t[k]/mi_t[0]<<endl;
    }


}

void sc_dos_tabulate1D(string fname,double* (*f)(double,double,double),double g0,double tp,int N)
{
    fstream outfile(path(fname,7),fstream::out);
    int npi=100000;
    double omegamin=0.;
    double omegamax=0.;
    double x=-pi;
    double step=2.0*pi/npi;

    for(int i=0;i<npi;i++)
    {
      double omega=ek(x,tp);
      double gamma=g(x,g0);
      double omp=omega+gamma;
      double omm=omega-gamma;
      if (omm < omegamin) omegamin=omm;
      if (omp > omegamax) omegamax=omp;
      x=x+step;
    }

    double h=(omegamax-omegamin)/(N-1);
    double* tmp;

    if(!outfile.good())
    {
        cout<<"nie otwarty plik!"<<endl;
    }

    double omega=omegamin;
    while(omega<omegamax)
    {
        omega=omega+h;
        tmp=f(omega,g0,tp);

        outfile<<std::setprecision(12)<<omega<<" "<<std::setprecision(12)<<tmp[0]<<" "<<std::setprecision(12)<<tmp[1]<<" "<<std::setprecision(12)<<tmp[2]<<endl;


    }

    outfile.close();
}
