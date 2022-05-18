/*
 *  Vectores.h
 *  
 *
 *  Created by Paulo Amorim on 14/Jul/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef VECTORES_H
#define VECTORES_H

#include <cmath>
#include <complex>
#include "Dados.h"
#include "Tridiagonal.h"

using namespace std;

typedef complex<double> Complex;

class vect {
public:
	Complex *comps;  //vector.comps[i] é a componente i
	int dim;
	
//	vect();
	vect(Dados *);		//Constr. com Dados!
	vect(int);
	vect(int , Complex);//args(int dimension) assign value to all the vector entries
	vect(const vect&);	//Copy constructor. Tem de ser!!
	~vect();
	
	void show() const;		//Mostrar vector
	
	vect &operator=(const  vect &);  //assignment 
	vect &operator+=(const vect &);  //
	vect &operator-=(const vect &);  //
	vect operator-() const;		//Unary minus
	friend vect operator+(const vect &, const vect &);
	friend vect operator+(const vect &,const Complex );//Vector plus a constant
	friend vect operator-(const vect &,const Complex );//Vector minus a constant 
	friend vect operator-(const vect &, const vect &);
	friend vect operator*(const Complex , const vect &);//Product BY CONSTANT
	friend Complex operator*(const vect &, const vect &);	//Scalar Product
	friend vect operator*(const vect &, const tridiagonal &); //Product of
																//Tridiagonal by vector
	
	vect realvect (void);
	vect conjugate (void);
	friend double MaxDif (const vect&, const vect&);
	friend vect Compare (const vect&, const vect&);	
	friend double NormaL2 (const vect&, double);

	//friend vect SolveTri(const tridiagonal &mat, const vect &uu);
	vect SolveTri(const tridiagonal& mat);
};





#endif VECTORES_H