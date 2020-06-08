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
    FILE *r= fopen("W.txt", "r");
    for (int i = 0; i < N-1; ++i)
    {
        fscanf(r, "%lf   %lf   %lf   %lf   %lf", &aux1,&aux2,&aux3,&aux4,&aux5);
        gamma_med=gamma_med+aux2;
    }
    fclose(r);
    return(gamma_med/N);
}

double trueKmed(int N)
{
    double aux1,aux2,aux3,aux4,aux5;
    double gamma_med=0;
    FILE *r= fopen("W.txt", "r");
    for (int i = 0; i < N-1; ++i)
    {
        fscanf(r, "%lf   %lf   %lf   %lf   %lf", &aux1,&aux2,&aux3,&aux4,&aux5);
        gamma_med=gamma_med+aux1;
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
    double F; 
    double K_1; 
    double gamma_1; 
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

    if(regime==2)
    {
        printf("F: ");
        std::cin >>F;
        printf("K_1: ");
        std::cin >>K_1;
        printf("gamma_1: ");
        std::cin >>gamma_1;
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
    double A_med;
    double w_prime;
    double c;
    double a1;

    int N;
    printf("N: ");
    std::cin >>N;
    if(regime==1)
    {
    	A_med_squared=pow((double)K_med/(N),2)/(1.0+pow(gamma_med,2));
    	printf("A_med=%lf\n",sqrt(A_med_squared) );
    } 
    if(regime==2)
    {
        c=(K_med/500.0)*1.0/sqrt(pow(K_med/500.0-1,2)+pow(gamma_med,2));
        printf("A_med/A_1=%lf\n",c);
        a1=(pow(1.0/F,2))*(pow(K_med*0.998-1-c*c*K_1*(1.0/K_med)*499*(K_med/500.0-1),2)+pow(gamma_1+gamma_med*499*K_1*(1.0/K_med)*c*c,2));
        a1=sqrt(1.0/a1);
        printf("A_1=%lf\n",a1 );
        A_med=a1*c;
        printf("A_med=%lf\n",A_med);
        w_prime=0.5*a1*a1*((gamma_1+499*K_1*(1.0/K_med)*c*c*gamma_med)/(gamma_1+499*K_1*(1.0/K_med)*gamma_med));
        printf("W'=%lf\n",w_prime );
    } 

	FILE *gplotpipe = popen( "gnuplot -persist", "w" );
    if(regime==0)
    {
    	if(W==0)
    	{
    		fprintf(gplotpipe, "Z(x,y)=%lf\n",gamma_med);
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
    if(regime==2)
    {
        if(W==1)
        {
            fprintf(gplotpipe, "Z(x,y)=%.15lf*0.5*%lf*0.998*(1/((0.998*y-1)*(0.998*y-1)+x*x))*(",A_med*A_med,K_med);
            fprintf(gplotpipe, "(x-%lf)*(0.998*y-1)-(y/%lf-1)*(x+0.998*y*%lf-%lf)",gamma_med,K_med,gamma_med,gamma_med);
            fprintf(gplotpipe, ")-%.15lf*%.15lf*%lf*%lf*(1.0/y)*(y/%lf-1-x/%lf+1)\n",w_prime,w_prime,gamma_med,K_med,K_med,gamma_med);
        }
        if(W==0)
        {
            fprintf(gplotpipe, "Z(x,y)=%.15lf*0.5*(",A_med*A_med);
            fprintf(gplotpipe, "%.15lf+(1/((0.998*y-1)*(0.998*y-1)+x*x))*((y/%lf-1)*x-(x-%lf)*(0.998*y-1)+x*%lf*%lf*(y/%lf-x/%lf))",gamma_med,K_med,gamma_med,gamma_med,gamma_med,K_med,gamma_med);
            fprintf(gplotpipe, ")+%lf*%.15lf*%.15lf*(x+%lf*(y/%lf-x/%lf))/y\n",K_med,w_prime,w_prime,gamma_med,K_med,gamma_med);
        }
        if(W==2)
        {
            fprintf(gplotpipe, "Z(x,y)=-%.15lf*0.5*x*%lf*(1/y)*(",A_med*A_med,K_med);
            fprintf(gplotpipe, "1+(1/((0.998*y-1)*(0.998*y-1)+x*x))*((y/%lf-1)*(y/%lf-1)+(%lf*(y/%lf-x/%lf))*(%lf*(y/%lf-x/%lf)))+2*(1/((0.998*y-1)*(0.998*y-1)+x*x))*(x*%lf*(y/%lf-x/%lf)-(y/%lf-1)*(%lf*0.998-1))",K_med,K_med,gamma_med,K_med,gamma_med,gamma_med,K_med,gamma_med,gamma_med,K_med,gamma_med,K_med,K_med);
            fprintf(gplotpipe, ")-x*%lf*(1/y)*%.15lf*%.15lf\n",K_med,w_prime,w_prime);
        }
    }

    //fprintf(gplotpipe, "set pm3d\n");
    fprintf(gplotpipe, "set surface\n");
    fprintf(gplotpipe, "set samples 200,200\n");
    fprintf(gplotpipe, "set isosamples 30,30\n");
    if(test ==1)
    {
        if(regime==1)
        {
            fprintf(gplotpipe, "Z1(x,y)=%lf*0.5*(",A_med_squared);
            fprintf(gplotpipe, "%lf-(x-%lf)*(x*%lf+0.998*y-1)/((0.998*y-1)*(0.998*y-1)+x*x)+(y/%lf-1)*x*(1+%lf*%lf)/((0.998*y-1)*(0.998*y-1)+x*x))\n",gamma_med,gamma_med,gamma_med,K_med,gamma_med,gamma_med);
            fprintf(gplotpipe, "Z2(x,y)=%lf*0.5*(",A_med_squared);
            fprintf(gplotpipe, "(x-%lf)*(0.998*y-1)*%lf*0.998/((0.998*y-1)*(0.998*y-1)+x*x)-0.998*(y-%lf)*(%lf*(0.998*y-1)+x)/((0.998*y-1)*(0.998*y-1)+x*x)",gamma_med,K_med,K_med,gamma_med);
            fprintf(gplotpipe, ")\n");
            fprintf(gplotpipe, "Z3(x,y)=%lf*0.5*(",A_med_squared);
            fprintf(gplotpipe, "-x*%lf/y+x*(x-%lf)*(x+%lf)*(%lf/y)/((0.998*y-1)*(0.998*y-1)+x*x)-(y/%lf-1)*(y/%lf-1)*(1+%lf*%lf)*x*%lf/(y*((0.998*y-1)*(0.998*y-1)+x*x))-2*x*(1-%lf/y)*(x*%lf-0.998*y+1)/((0.998*y-1)*(0.998*y-1)+x*x)",K_med,gamma_med,gamma_med,K_med,K_med,K_med,gamma_med,gamma_med,K_med,K_med,gamma_med);
            fprintf(gplotpipe, ")\n");
            fprintf(gplotpipe, "splot [%lf:%lf] [%lf:%lf] (Z1(x,y)+Z2(x,y)+Z3(x,y))*%lf, 'W.txt' u 1:2:($3+$4+$5)\n",gammarangemin,gammarangemax,Krangemin,Krangemax,cte);
        }
        if(regime==2)
        {
            fprintf(gplotpipe, "Z1(x,y)=%.15lf*0.5*%lf*0.998*(1/((0.998*y-1)*(0.998*y-1)+x*x))*(",A_med*A_med,K_med);
            fprintf(gplotpipe, "(x-%lf)*(0.998*y-1)-(y/%lf-1)*(x+0.998*y*%lf-%lf)",gamma_med,K_med,gamma_med,gamma_med);
            fprintf(gplotpipe, ")-%.15lf*%.15lf*%lf*%lf*(1.0/y)*(y/%lf-1-x/%lf+1)\n",w_prime,w_prime,gamma_med,K_med,K_med,gamma_med);
            fprintf(gplotpipe, "Z2(x,y)=%.15lf*0.5*(",A_med*A_med);
            fprintf(gplotpipe, "%.15lf+(1/((0.998*y-1)*(0.998*y-1)+x*x))*((y/%lf-1)*x-(x-%lf)*(0.998*y-1)+x*%lf*%lf*(y/%lf-x/%lf))",gamma_med,K_med,gamma_med,gamma_med,gamma_med,K_med,gamma_med);
            fprintf(gplotpipe, ")+%lf*%.15lf*%.15lf*(x+%lf*(y/%lf-x/%lf))/y\n",K_med,w_prime,w_prime,gamma_med,K_med,gamma_med);
            fprintf(gplotpipe, "Z3(x,y)=-%.15lf*0.5*x*%lf*(1/y)*(",A_med*A_med,K_med);
            fprintf(gplotpipe, "1+(1/((0.998*y-1)*(0.998*y-1)+x*x))*((y/%lf-1)*(y/%lf-1)+(%lf*(y/%lf-x/%lf))*(%lf*(y/%lf-x/%lf)))+2*(1/((0.998*y-1)*(0.998*y-1)+x*x))*(x*%lf*(y/%lf-x/%lf)-(y/%lf-1)*(%lf*0.998-1))",K_med,K_med,gamma_med,K_med,gamma_med,gamma_med,K_med,gamma_med,gamma_med,K_med,gamma_med,K_med,K_med);
            fprintf(gplotpipe, ")-x*%lf*(1/y)*%.15lf*%.15lf\n",K_med,w_prime,w_prime);
            fprintf(gplotpipe, "splot [%lf:%lf] [%lf:%lf] (Z1(x,y)+Z2(x,y)+Z3(x,y))*%lf, 'W.txt' u 1:2:($3+$4+$5)\n",gammarangemin,gammarangemax,Krangemin,Krangemax,cte);
        }
    }
    if(test!=1)
    {
        fprintf(gplotpipe, "splot [%lf:%lf] [%lf:%lf] Z(x,y)*%lf w l palette, 'W.txt' u 2:1:%d w p palette pt 7 ps 1\n",gammarangemin,gammarangemax,Krangemin,Krangemax,cte,W+3);
    }
    printf("<g>=%lf   <K>=%lf\n",truegammamed(N),trueKmed(N));

	return 0;
}