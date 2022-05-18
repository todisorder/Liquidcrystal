/*
 *  Tridiagonal.h
 *  
 *
 *  Created by Paulo Amorim on 13/Jul/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef TRIDIAGONAL_H
#define TRIDIAGONAL_H

#include <cmath>
#include "Dados.h"
//#include "Vectores.h"
#include <complex>

using namespace std;

typedef complex<double> Complex;
//typedef double Complex;


class tridiagonal
{
public:
	
	Complex *ld, *dd, * ud;
	int dim;

	tridiagonal(int);			//Construtor: tridiagonal T (4); deve criar mat. tridiagonal T 4x4.
	tridiagonal(Dados *);	//Const. com Dados
	~tridiagonal();
	
	void set_lowdiag(Complex *);
	void set_diagdiag(Complex *);
	void set_updiag(Complex *);
	
	void show() const;			//Mostrar matriz
	
	tridiagonal &operator=(const tridiagonal &);  //assignment tri1 = tri2
	tridiagonal &operator+=(const tridiagonal &);  //
	tridiagonal &operator-=(const tridiagonal &);  //
	tridiagonal operator-() const;		//Unary minus
	friend tridiagonal operator+(const tridiagonal &, const tridiagonal &);
	friend tridiagonal operator-(const tridiagonal &, const tridiagonal &);
	friend tridiagonal operator*(const Complex , const tridiagonal &);//Product BY CONSTANT

};	
	
#endif   //TRIDIAGONAL_H