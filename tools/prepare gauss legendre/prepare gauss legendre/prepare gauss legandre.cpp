

// prepare gauss legandre.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cstdlib>
#include <complex>
#include <Eigenvalues>
#include <Dense>
#include <iomanip>
#include <sstream>

struct Polynomial
{
    unsigned int degree;
    double* coef;
};

Polynomial zero_init(unsigned int degree)
{
    Polynomial res = { degree , 0 };
    res.coef = new double[res.degree+1];

    for (int i = 0; i <= res.degree; i++)
        res.coef[i] = 0.0;

    return res;
}

Polynomial mul(Polynomial a, Polynomial b)
{
    Polynomial res = { a.degree + b.degree, 0 };
    res.coef = new double[res.degree+1];

    for (int i = 0; i <= res.degree; i++)
        res.coef[i] = 0.0;

    for (int i = 0; i <= a.degree; i++)
    {
        for (int j = 0; j <= b.degree; j++)
        {
            res.coef[i + j] += a.coef[i] * b.coef[j];
        }
    }

    return res;
}

Polynomial add(Polynomial a, Polynomial b)
{
    Polynomial res = { std::max(a.degree,b.degree), 0 };
    res.coef = new double[res.degree+1];

    for (int i = 0; i <= std::min(a.degree,b.degree); i++)
        res.coef[i] = a.coef[i]+b.coef[i];

    for (int i = std::min(a.degree, b.degree)+1; i <= res.degree; i++)
    {
        if (a.degree > b.degree) res.coef[i] = a.coef[i];
        else res.coef[i] = b.coef[i];
    }

    return res;
}

Polynomial df(Polynomial a)
{
    Polynomial res = zero_init(a.degree - 1);

    for (int i = 0; i <= res.degree; i++)
    {
        res.coef[i] = (i + 1) * a.coef[i + 1];
    }

    return res;
}

double eval(Polynomial a, double x)
{
    double res = 0.0;
    double xi = 1.0;
    for (int i = 0; i <= a.degree; i++)
    {
        res += xi * a.coef[i];
        xi *= x;
    }

    return res;

}
void print(Polynomial a)
{
    bool first = true;
    for (int i = 0; i <= a.degree; i++)
    {
        if (a.coef[i] != 0.0)
        {
            if (!first) std::cout << " + ";
            else first = false;

            std::cout << a.coef[i];
            if (i == 1) std::cout << "x";
            if (i >  1) std::cout<< "x^" << i;
        }
    }
    std::cout << std::endl;
}

Polynomial Legendre(int n)
{
    Polynomial Lpp = zero_init(n);
    Lpp.coef[0] = 1.0;

    if (n == 0) return Lpp;

    Polynomial Lp = zero_init(n);
    Lp.coef[1] = 1.0;

    if (n == 1) return Lp;
    
    Polynomial help1 = zero_init(1);
    Polynomial help0 = zero_init(0);
    Polynomial res = zero_init(n);

    for (int i = 2; i <= n; i++)
    {
        help1.coef[1] = (2.0 * i - 1.0) / double(i);
        help0.coef[0] = - (i - 1.0) / double(i);

        res = add(mul(help1, Lp), mul(help0, Lpp));
        Lpp = Lp;
        Lp = res;
    }

    return res;
}

void prepare_gl(int n, double ak[], double xk [])
{
    double y;
    Eigen::MatrixXd compmat(n,n);
    Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces;
    std::complex <double> temp(0.0,0.0);
   
    Polynomial L = Legendre(n);
    Polynomial L1 = Legendre(n+1);

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            compmat(i,j)=0;
    
    for(int i=0;i<n-1;i++)
    {
        compmat(i+1,i)=1.0;
        compmat(i,n-1)=-L.coef[i]/L.coef[n];
    }

    ces.compute(compmat);

    Polynomial dL = df(L);

    for(int i=0;i<n;i++)
    {
        temp=ces.eigenvalues()[i];
        xk[i]=temp.real();
       // std::cout <<"xk["<<i<<"] = "<< xk[i] << " L(xk["<<i<<"]) = " << eval(L, xk[i]) << std::endl;

        y=(n+1)*eval(L1,xk[i]) * eval(dL,xk[i]);
        if(abs(y)>1.e-10)
            ak[i]=-2.0/y;
    }

    free(dL.coef);
    free(L1.coef);
    free(L.coef);
}

int main()
{
    double* ak;
    double* xk;

    std::cout << "//gl.h" << std::endl;
    std::cout << "#ifndef GL_CONSTANTS" << std::endl;
    std::cout << "#define GL_CONSTANTS" << std::endl;

    int k = 10;
    std::cout << "const int k =" << k << ";" << std::endl;
    std::cout << "__constant__ int d_k =" << k << ";" << std::endl;
    ak = new double[k];
    xk = new double[k];
    prepare_gl(k, ak, xk);
    
    std::cout << "__constant__ double ak[] = {";
    for (int i = 0; i < k; i++)
    {
        std::cout << std::setprecision(16) << ak[i];
        if(i !=k-1) std::cout<< ", ";
    }
    std::cout << "};" << std::endl;

    std::cout << "__constant__ double xk[] = {";
    for (int i = 0; i < k; i++)
    {
        std::cout << std::setprecision(16) << xk[i];
        if(i !=k-1) std::cout<< ", ";

    }
    std::cout << "};" << std::endl;

    std::cout << "#endif" << std::endl;
    free(ak);
    free(xk);
    return 0;
   
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
