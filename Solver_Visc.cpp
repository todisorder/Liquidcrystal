/*
 *  Solver_Visc.cpp
 *  
 *
 *  Created by Paulo Amorim on 2/Aug/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "Solver_Visc.h"
#include <iostream>

using namespace std;

/********************************************************************/
//			Construtor com data
/********************************************************************/
Solver_Visc::Solver_Visc(Dados *data)
: UU(data), RR(data), WW(data), WWini(data), XX(data), SOL_EXACTAUU(data)
{		//Viste este sintaxe para inicializar? está em Capper, pág. 292.
		//é essencial ter um construtor com (data) para estes UU, etc!

		
	
	int nn = data->get_int("numpt");
	int numiter = data->get_int("numiter");
	double dx = data->get_double("dx");
	double dt = data->get_double("dt");
	double xmin = data->get_double("xmin");
	double focus = data->get_double("focus");
	double eps = data->get_double("eps");
	double delta = data->get_double("delta");
	numpt = nn;

	double tfinal = dt*numiter;


	/********************************************************************/
	//			Vector espacial XX
	/********************************************************************/
	for (int i=0; i<numpt; i++) {
		XX.comps[i] = xmin + i*dx;
	}
	Complex xmid = 0.5*(XX.comps[0] + XX.comps[numpt-1]);		//Pto médio do intervalo
	/********************************************************************/
	//			Dado inicial para UU -- Tipo solitão
	/********************************************************************/
	
	
	double a=.5, c=5., d = -0., b=1.;		//Parametros para Solitão
	double a1=0., c1=3., d1 = -10., b1=1.;		//Parametros para Solitão
	//
	//double eps = .05;
	double lmda = 1. , cc = 3./2.  , EE = lmda - 0.25*cc*cc , gam = 1.,	ac = 1./(gam - cc);	
	Complex CONST = sqrt(((2.*EE)/abs(ac+focus)));
	for (int j = 0; j <= numpt-1; j++)
    {
/*
			UU.comps[j] = eps*sqrt(1.5)*(1./(cosh(0.5*eps*sqrt(3.)*XX.comps[j])))
			*exp(Complex(0.,1.)*(0.5*eps+1.)*XX.comps[j]);
		UU.comps[j] = abs(UU.comps[j]);
*/
    }
	
	for (int j = 0; j <= numpt-1; j++)
	{
//		Solucao que tende para zero em tempo finito
		//		UU.comps[j] = 0.5*XX.comps[j]*(1.-XX.comps[j]);
	}

	for (int j = 0; j <= numpt-1; j++)
	{
		//		Solucao que tende para chapéu em tempo finito
//		UU.comps[j] = min( abs(0.5*XX.comps[j]*XX.comps[j]),abs(1.-XX.comps[j])) ;
		// TIPO solitão:
//        UU.comps[j] = 1.*norm(CONST*(1./(cosh(sqrt(100.)*(XX.comps[j] - 1.)))));
        UU.comps[j] = 1.* exp(Complex(0.,1.)*20.*XX.comps[j])*norm(1.*(1./(cosh(sqrt(5.)*(XX.comps[j] + 1.)))));
	}
	
	//Sol "exacta" q_A
	//O SEGUINTE NAO INTERESSA.
	for (int j = 0; j <= numpt-1; j++){
		SOL_EXACTAUU.comps[j] = eps*exp(Complex(0.,1.)*((0.5*eps+1.)*XX.comps[j]+tfinal*(1.-eps*eps+eps)))*sqrt(1.5)*(1./(cosh(0.5*eps*sqrt(3.)*(XX.comps[j]+2.*tfinal))));
	}
	
	/********************************************************************/
	//			Solução exacta de UU (Quando há)
	/********************************************************************/
	for (int j = 0; j <= numpt-1; j++)
	{		
		//********* Sol. exacta para Laurençot (fluxo linear)
	//	SOL_EXACTAUU.comps[j] = CONST*(1./(cosh(sqrt(EE)*(XX.comps[j] - cc*tfinal ))))
	//	*exp(Complex(0.,1.)*0.5*cc*(XX.comps[j] - cc*tfinal ))*exp(Complex(0.,1.)*lmda*Complex(tfinal));
	}
	
	/********************************************************************/
	//			Dado inicial para RR -- Tipo solitão
	/********************************************************************/
	
	for (int j = 0; j <= numpt-1; j++)
	{
		// TIPO solitão:
		RR.comps[j] = 1.*norm(1.*(1./(cosh(sqrt(10.)*(XX.comps[j] - 0.)))));
	}
	
	
	/********************************************************************/
	//			Dado inicial para WW -- Tipo solitão
	/********************************************************************/
	
	
	for (int j = 0; j <= numpt-1; j++)
	{
		// TIPO solitão:
		WW.comps[j] = .3*norm(1.*(1./(cosh(sqrt(100.)*(XX.comps[j] - 0.)))));
//        WWini.comps[j] = WW.comps[j];
	}
	
	
	//Dado inicial tipo BOX
	//Runge - Kutta não dá para isto!!!!!!
/*
	double aa=10., dd = -1., bb=1.;
	int indice1 = 15*(numpt-1)/30;
	int indice2 = 20*(numpt-1)/30;
	Complex xter = XX.comps[indice1];
	Complex x2ter = XX.comps[indice2];
	for (int j=0; j<=indice1; j++) 
	{
			RR.comps[j] = bb * exp(-100.*(XX.comps[j] - xter)*(XX.comps[j] - xter));
			WW.comps[j] = bb * exp(-100.*(XX.comps[j] - xter)*(XX.comps[j] - xter));
    }	
	for (int j=indice1+1; j<=indice2; j++) {
			RR.comps[j] = Complex(bb,0.);
			WW.comps[j] = Complex(bb,0.);
	}
	for (int j=indice2+1; j<=numpt-1; j++)
    {
			RR.comps[j] = bb * exp(-100.*(XX.comps[j] - x2ter)*(XX.comps[j] - x2ter));
			WW.comps[j] = bb * exp(-100.*(XX.comps[j] - x2ter)*(XX.comps[j] - x2ter));
    }
*/
	
	
	
}

/********************************************************************/
//				Método de RK 4 Explícito
//
/********************************************************************/

Complex Flux(const Complex v, const double a, const double b, const double c)
{
    double lam = (2/3)*c*(a-b);
	Complex tempo = a*v + lam*v*v*v;
	return tempo;		
}

///********************************************************************/
Complex FuncRKUU(Dados *data, const Complex ujm, const Complex uj, const Complex ujj, const Complex rj, const Complex xxj, double tcurrent)  // 6 --> 7
{
	double dx = data->get_double("dx");
    double aa = data->get_double("aa");
    double HH = data->get_double("HH");

	Complex iii = Complex(0.,1.);
	Complex tempo = iii*( (1./(dx*dx))*(ujm - 2.*uj + ujj) + Complex(1.,0.)*rj*uj - aa*uj*norm(uj) - HH*HH * xxj*xxj*uj  );
	
	return 1.*tempo;	
}
///********************************************************************/
Complex FuncRKRR(Dados *data, const Complex wj, double tcurrent)    // 4 --> 3
{
	Complex tempo = wj;
	return tempo;	
}
///********************************************************************/
Complex FuncRKWW(Dados *data, const Complex rjm, const Complex rj, const Complex rjj, const Complex wjm, const Complex wj, const  Complex wjj, const  Complex uj, double tcurrent)  // +1 -2 -1= -2
{
	double dx = data->get_double("dx");
	double delta = data->get_double("delta");	
	double eps = data->get_double("eps");
    double bb = data->get_double("bb");
    double alpha = data->get_double("alpha");
    double beta = data->get_double("beta");
    double gamma = data->get_double("gamma");
    
    Complex Fmais = Flux((rjj-rj)/dx,alpha,beta,gamma);
    Complex Fmenos = Flux((rj-rjm)/dx,alpha,beta,gamma);
    
	Complex tempo = Complex(1.,0.)*(1./dx)*(Fmais - Fmenos) +  Complex(1.,0.)*(eps/(dx*dx))*(wjm - 2.*wj + wjj) - Complex(1.,0.)*bb*rj + Complex(1.,0.)*norm(uj);
	
	return Complex(1.,0.)*tempo;
}

void Solver_Visc::RKVISC(Dados *data, double tcurrent){
	
	/********************************************************************/
	//		Buscar dados
	/********************************************************************/
	int dim = data->get_int("numpt");
	double dt = data->get_double("dt");
	double dx = data->get_double("dx");

	vect unk (UU);
	vect un (UU);
	vect kapaum(dim);
	vect kapadois(dim);
	vect kapatres(dim);
	vect kapaquatro(dim);
	
    int j = 1;
    vect FF(dim);
	
	/********************************************************************/
	/*******************************		UU		*********************/
	for (int i=1; i<dim-2; i++) {
		kapaum.comps[i] = dt*(FuncRKUU(data, UU.comps[i-1],UU.comps[i],UU.comps[i+1],RR.comps[i],XX.comps[i],tcurrent));
	}
	
	for (int i=1; i<dim-2; i++) {
		kapadois.comps[i] = dt*(FuncRKUU(data, UU.comps[i-1]+0.5*kapaum.comps[i-1],UU.comps[i]+0.5*kapaum.comps[i],UU.comps[i+1]+0.5*kapaum.comps[i+1],RR.comps[i],XX.comps[i],tcurrent));
	}
	
	for (int i=1; i<dim-2; i++) {
		kapatres.comps[i] = dt*(FuncRKUU(data, UU.comps[i-1]+0.5*kapadois.comps[i-1],UU.comps[i]+0.5*kapadois.comps[i],UU.comps[i+1]+0.5*kapadois.comps[i+1],RR.comps[i],XX.comps[i],tcurrent));
	}
	
	for (int i=1; i<dim-2; i++) {
		kapaquatro.comps[i] = dt*(FuncRKUU(data, UU.comps[i-1]+kapatres.comps[i-1],UU.comps[i]+kapatres.comps[i],UU.comps[i+1]+kapatres.comps[i+1],RR.comps[i],XX.comps[i],tcurrent));
	}
	
	for (int i=1; i<dim-2; i++) {
		UU.comps[i] = UU.comps[i] + (1./6.)*(kapaum.comps[i] + 2.*kapadois.comps[i] + 2.*kapatres.comps[i] + kapaquatro.comps[i]);
	}
	
	/********************************************************************/
	/************************		RR		*****************************/
	//RK4 reduz-se a Euler se F não depende de RR...	//
	for (int i=1; i<dim-2; i++) {
		RR.comps[i] = RR.comps[i] + dt*FuncRKRR(data, WW.comps[i], tcurrent);
	}
	/********************************************************************/
	/********************************		WW		*********************/
	for (int i=1; i<dim-2; i++) {
		kapaum.comps[i] = dt*(FuncRKWW(data, RR.comps[i-1], RR.comps[i], RR.comps[i+1], WW.comps[i-1], WW.comps[i], WW.comps[i+1], UU.comps[i], tcurrent));
	}
	
	for (int i=1; i<dim-2; i++) {
		kapadois.comps[i] = dt*(FuncRKWW(data, RR.comps[i-1], RR.comps[i], RR.comps[i+1], WW.comps[i-1]+0.5*kapaum.comps[i-1], WW.comps[i]+0.5*kapaum.comps[i], WW.comps[i+1]+0.5*kapaum.comps[i+1], UU.comps[i], tcurrent));
	}
	
	for (int i=1; i<dim-2; i++) {
		kapatres.comps[i] = dt*(FuncRKWW(data, RR.comps[i-1], RR.comps[i], RR.comps[i+1], WW.comps[i-1]+0.5*kapadois.comps[i-1], WW.comps[i]+0.5*kapadois.comps[i], WW.comps[i+1]+0.5*kapadois.comps[i+1], UU.comps[i], tcurrent));
	}
	
	for (int i=1; i<dim-2; i++) {
		kapaquatro.comps[i] = dt*(FuncRKWW(data, RR.comps[i-1], RR.comps[i], RR.comps[i+1], WW.comps[i-1]+kapatres.comps[i-1], WW.comps[i]+kapatres.comps[i], WW.comps[i+1]+kapatres.comps[i+1], UU.comps[i], tcurrent));
	}
	
	for (int i=1; i<dim-2; i++) {
		WW.comps[i] = WW.comps[i] + (1./6.)*(kapaum.comps[i] + 2.*kapadois.comps[i] + 2.*kapatres.comps[i] + kapaquatro.comps[i]);
	}
	
}


/********************************************************************/
//		Resolver o problema Schrodinger Linear (para testar)
/********************************************************************/
//void Solver_Visc::Linear(Dados *data)
//{
//    int dim = data->get_int("numpt");
//    vect un (UU*MAT_RHS);
//    un.SolveTri(MAT);  //deve bastar isto... muda un.
//    UU = un;
//}

