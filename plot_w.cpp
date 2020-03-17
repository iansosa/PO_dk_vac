#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>
#include <cassert>
#include <vector>
#include <iostream>
#include <armadillo>
#include <utility>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

double truegammamed(int N)
{
    double aux1,aux2,aux3,aux4,aux5;
    double gamma_med=0;
    FILE *r= fopen("TempProm.txt", "r");
    for (int i = 0; i < N-1; ++i)
    {
        fscanf(r, "%lf   %lf   %lf   %lf   %lf", &aux1,&aux2,&aux3,&aux4,&aux5);
        gamma_med=gamma_med+aux1;
    }
    fclose(r);
    return(gamma_med/N);
}

double trueKmed(int N)
{
    double aux1,aux2,aux3,aux4,aux5;
    double gamma_med=0;
    FILE *r= fopen("TempProm.txt", "r");
    for (int i = 0; i < N-1; ++i)
    {
        fscanf(r, "%lf   %lf   %lf   %lf   %lf", &aux1,&aux2,&aux3,&aux4,&aux5);
        gamma_med=gamma_med+aux2;
    }
    fclose(r);
    return(gamma_med/N);
}

int main()
{
	using namespace boost::numeric::odeint;

    int test;
    printf("Test? (1 YES)");
    std::cin >>test; 

    int regime;
    printf("Regime: [1a(0), 1b(1), 2(2)]: ");
    std::cin >>regime; 
    if(regime==0)
    {
    	printf("1a\n");
    }
    if(regime==1)
    {
    	printf("1b\n");
    }
    if(regime==2)
    {
    	printf("2\n");
    }

    int W;
    if(test != 1)
    {
        printf("Transfer: [1(0), omega(1), gamma(2)]: ");
        std::cin >>W; 
        if(W==0)
        {
            printf("W1\n");
        }
        if(W==1)
        {
            printf("WOmega\n");
        }
        if(W==2)
        {
            printf("WGamma\n");
        }
    }

    double gamma_med;
    printf("<gamma>: ");
    std::cin >>gamma_med; 

    double K_med;
    printf("<K>: ");
    std::cin >>K_med;

    double gammarangemin;
    double gammarangemax;
    printf("G[x,]: ");
    std::cin >>gammarangemin; 
    printf("G[,x]: ");
    std::cin >>gammarangemax; 
    double Krangemin;
    double Krangemax;
    printf("K[x,]: ");
    std::cin >>Krangemin; 
    printf("K[,x]: ");
    std::cin >>Krangemax; 

    double cte;
    printf("MC: ");
    std::cin >>cte;


    double A_med_squared;
    int N;
    if(regime==1)
    {
    	printf("N: ");
    	std::cin >>N;

    	A_med_squared=pow((double)K_med/(N),2)/(1.0+pow(gamma_med,2));
    	printf("A_med=%lf\n",sqrt(A_med_squared) );
    } 

	FILE *gplotpipe = popen( "gnuplot -persist", "w" );
    if(regime==0)
    {
    	if(W==0)
    	{
    		fprintf(gplotpipe, "Z(x,y)=-%lf\n",gamma_med);
    	}
    	if(W==1)
    	{
    		fprintf(gplotpipe, "Z(x,y)=-(1-%lf/y)*%lf",K_med,gamma_med);
    		fprintf(gplotpipe, "+(x-%lf)*%lf/y\n",gamma_med,K_med);
    	}
    	if(W==2)
    	{
    		fprintf(gplotpipe, "Z(x,y)=-x*%lf/y\n",K_med);
    	}
    }
    if(regime==1)
    {
    	if(W==0)
    	{
    		fprintf(gplotpipe, "Z(x,y)=%lf*0.5*(",A_med_squared);
    		fprintf(gplotpipe, "%lf-(x-%lf)*(x*%lf+0.998*y-1)/((0.998*y-1)*(0.998*y-1)+x*x)+(y/%lf-1)*x*(1+%lf*%lf)/((0.998*y-1)*(0.998*y-1)+x*x))\n",gamma_med,gamma_med,gamma_med,K_med,gamma_med,gamma_med);
    	}
    	if(W==1)
    	{
    		fprintf(gplotpipe, "Z(x,y)=%lf*0.5*(",A_med_squared);
    		fprintf(gplotpipe, "(x-%lf)*(0.998*y-1)*%lf*0.998/((0.998*y-1)*(0.998*y-1)+x*x)-0.998*(y-%lf)*(%lf*(0.998*y-1)+x)/((0.998*y-1)*(0.998*y-1)+x*x)",gamma_med,K_med,K_med,gamma_med);
    		fprintf(gplotpipe, ")\n");
    	}
    	if(W==2)
    	{
    		fprintf(gplotpipe, "Z(x,y)=%lf*0.5*(",A_med_squared);
    		fprintf(gplotpipe, "-x*%lf/y+x*(x-%lf)*(x+%lf)*(%lf/y)/((0.998*y-1)*(0.998*y-1)+x*x)-(y/%lf-1)*(y/%lf-1)*(1+%lf*%lf)*x*%lf/(y*((0.998*y-1)*(0.998*y-1)+x*x))-2*x*(1-%lf/y)*(x*%lf-0.998*y+1)/((0.998*y-1)*(0.998*y-1)+x*x)",K_med,gamma_med,gamma_med,K_med,K_med,K_med,gamma_med,gamma_med,K_med,K_med,gamma_med);
    		fprintf(gplotpipe, ")\n");
    	}
    }

    fprintf(gplotpipe, "set pm3d\n");
    fprintf(gplotpipe, "set samples 200,200\n");
    fprintf(gplotpipe, "set isosamples 30,30\n");
    if(test ==1)
    {

            fprintf(gplotpipe, "Z1(x,y)=%lf*0.5*(",A_med_squared);
            fprintf(gplotpipe, "%lf-(x-%lf)*(x*%lf+0.998*y-1)/((0.998*y-1)*(0.998*y-1)+x*x)+(y/%lf-1)*x*(1+%lf*%lf)/((0.998*y-1)*(0.998*y-1)+x*x))\n",gamma_med,gamma_med,gamma_med,K_med,gamma_med,gamma_med);
            fprintf(gplotpipe, "Z2(x,y)=%lf*0.5*(",A_med_squared);
            fprintf(gplotpipe, "(x-%lf)*(0.998*y-1)*%lf*0.998/((0.998*y-1)*(0.998*y-1)+x*x)-0.998*(y-%lf)*(%lf*(0.998*y-1)+x)/((0.998*y-1)*(0.998*y-1)+x*x)",gamma_med,K_med,K_med,gamma_med);
            fprintf(gplotpipe, ")\n");
            fprintf(gplotpipe, "Z3(x,y)=%lf*0.5*(",A_med_squared);
            fprintf(gplotpipe, "-x*%lf/y+x*(x-%lf)*(x+%lf)*(%lf/y)/((0.998*y-1)*(0.998*y-1)+x*x)-(y/%lf-1)*(y/%lf-1)*(1+%lf*%lf)*x*%lf/(y*((0.998*y-1)*(0.998*y-1)+x*x))-2*x*(1-%lf/y)*(x*%lf-0.998*y+1)/((0.998*y-1)*(0.998*y-1)+x*x)",K_med,gamma_med,gamma_med,K_med,K_med,K_med,gamma_med,gamma_med,K_med,K_med,gamma_med);
            fprintf(gplotpipe, ")\n");
            fprintf(gplotpipe, "splot [%lf:%lf] [%lf:%lf] (Z1(x,y)+Z2(x,y)+Z3(x,y))*%lf, 'TempProm.txt' u 1:2:($3+$4+$5)\n",gammarangemin,gammarangemax,Krangemin,cte,Krangemax);
    }
    if(test!=1)
    {
        fprintf(gplotpipe, "splot [%lf:%lf] [%lf:%lf] Z(x,y)*%lf, 'TempProm.txt' u 1:2:%d\n",gammarangemin,gammarangemax,Krangemin,Krangemax,cte,W+3);
    }
    printf("<g>=%lf   <K>=%lf\n",truegammamed(N),trueKmed(N));

	return 0;
}