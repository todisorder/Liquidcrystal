/*
 *  Solver_Visc.h
 *  
 *
 *  Created by Paulo Amorim on 2/Aug/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SOLVER_VISC_H
#define SOLVER_VISC_H

#include <cmath>
#include <complex>
#include "Dados.h"
#include "Vectores.h"
using namespace std;

typedef complex<double> Complex;


class Solver_Visc
{
public:

	Solver_Visc(Dados *);		//Construtor 
	
	int numpt;

	vect UU;				//Sol de Schrodinger
	vect RR;				//
	vect WW;				//
	vect WWini;

	
	vect XX;				//Spatial Grid
//    tridiagonal MAT;        //Matriz A dos coeficientes do problema Schrodinger
//    tridiagonal MAT_NEWTON;    //Matriz DG para Newton DG(ukk - uk) = -G(uk), DG = A + DF
//    tridiagonal MAT_RHS;    //Matriz B, onde G = A + B + F
//	tridiagonal DF;			//Experiencia de memória
//    vect RHS_NEWTON;        //Vector -G do rhs de Newton
	vect SOL_EXACTAUU;

//    friend Complex FuncRKUU(Dados *, const Complex, const Complex, const Complex, const Complex, double tcurrent);
//    friend Complex FuncRKVV(Dados *, const Complex, const Complex, double tcurrent);
//    friend Complex FuncRKWW(Dados *, const Complex, const Complex, const Complex, const Complex, const Complex, const Complex, const Complex, const Complex, const Complex, double tcurrent);
//
//    friend Complex KK(const double, const Complex); // KK(eps, UU)
//    friend Complex GG(const    Complex);
//    friend Complex Flux(const Complex);
//    friend Complex Source(const Complex, double tcurrent);
	
	
	void RKVISC (Dados *, double tcurrent);
//    void Linear (Dados *);
		
};

#endif  //SOLVER_VISC_H
