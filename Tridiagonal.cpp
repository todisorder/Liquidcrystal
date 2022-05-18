/*
 *  Tridiagonal.cpp
 *  
 *
 *  Created by Paulo Amorim on 13/Jul/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>

#include "Tridiagonal.h"

using namespace std;

tridiagonal::tridiagonal(int d){		//Inicializar com zeros.
	dim = d;
	ld = new Complex  [dim - 1];
	dd = new Complex  [dim];
	ud = new Complex  [dim - 1];
	for (int j=0; j<dim-2; j++) {
		ld[j] = 0.0; ud[j] = 0.0; dd[j] = 0.0;		//Não sei se isto fica Complex!!!!!!
	}
	dd[dim-1] = 0.0;
}

tridiagonal::tridiagonal(Dados *data){		//Inicializar com Dados
	data->set();
	int nn = data->get_int("numpt");
	dim = nn;
	ld = new Complex  [dim - 1];
	dd = new Complex  [dim];
	ud = new Complex  [dim - 1];
	for (int j=0; j<dim-2; j++) {
		ld[j] = 0.0; ud[j] = 0.0; dd[j] = 0.0;		//Não sei se isto fica Complex!!!!!!
	}
	dd[dim-1] = 0.0;	
}




tridiagonal::~tridiagonal(){
	delete[] ld;
	delete[] dd;
	delete[] ud;	
//	cout << "Hello"<<endl;		//nunca aparece!!!!!
}



//------------------------------------------//------------------------------------------
//Build
//------------------------------------------//------------------------------------------
void tridiagonal::set_lowdiag(Complex *vect) {
	for (int j=0; j<=dim-2; j++) {	//d-2
		ld[j] = vect[j];
	}
}

void tridiagonal::set_diagdiag(Complex *vect) {
	for (int j=0; j<=dim-1; j++) {  //atenção, d-1.
		dd[j] = vect[j];
	}
}

void tridiagonal::set_updiag(Complex *vect) {
	for (int j=0; j<=dim-2; j++) {	//d-2
		ud[j] = vect[j];
	}
}

//------------------------------------------//------------------------------------------
//Show
//------------------------------------------//------------------------------------------
void tridiagonal::show() const {
	int aux = 0;
	Complex ZZ = 0.;
	cout << dd[0] <<"\t"<< ud[0]<<"\t";
	for (int i=2; i<=dim-1 ; i++) {
		cout << ZZ<<"\t";		//0th Line
	}
	cout <<endl;
	for (int line=1; line<=dim-2; line++) {	
		for (int col=0; col<=line-2; col++) {
			cout << ZZ <<"\t";
			++aux;
		}
		cout << ld[line-1] <<"\t" << dd[line]<<"\t" << ud[line]<<"\t";
		for (int col=aux+3; col<dim; col++) {
			aux =0.;
			cout << ZZ <<"\t";
		}
		cout << endl;
	}							//end of line 
	for (int i=0; i<dim-2; i++) {		//Last line:
		cout << ZZ <<"\t";
	}
	cout << ld[dim-2]<<"\t" << dd[dim-1]<<"\t" << endl;
}
	
//------------------------------------------//------------------------------------------
//Operator Overloads
//------------------------------------------//------------------------------------------

tridiagonal &tridiagonal::operator=(const tridiagonal &tri)
{
	dim = tri.dim;
	for (int i=0; i<=dim-2; i++) 
	{
		ld[i]=tri.ld[i]; dd[i]=tri.dd[i]; ud[i]=tri.ud[i];
	}
	dd[dim-1]=tri.dd[dim-1];
	
	
	//ld = tri.ld;//NÃO!!! assim se mudar um mudo o outro!!!!
	//dd = tri.dd;
	//ud = tri.ud;
	return *this;
}

tridiagonal &tridiagonal::operator+=(const tridiagonal &tri)
{
	dim = tri.dim;
	for (int i=0; i<=dim-2; i++) {
		ld[i]+=tri.ld[i]; dd[i]+=tri.dd[i]; ud[i]+=tri.ud[i];
	}
	dd[dim-1]+=tri.dd[dim-1];
	return *this;
}

tridiagonal &tridiagonal::operator-=(const tridiagonal &tri)
{
	dim = tri.dim;
	for (int i=0; i<=dim-2; i++) {
		ld[i]-=tri.ld[i]; dd[i]-=tri.dd[i]; ud[i]-=tri.ud[i];
	}
	dd[dim-1]-=tri.dd[dim-1];
	return *this;
}

tridiagonal tridiagonal::operator-() const
{
	tridiagonal tri(dim);
	for (int i=0; i<=dim-2; i++) {
		tri.ld[i]=-ld[i]; tri.dd[i]=-dd[i]; tri.ud[i]=-ud[i];
	}
	tri.dd[dim-1]=-dd[dim-1];
	return tri;
}

tridiagonal operator+(const tridiagonal &t1, const tridiagonal &t2)
{
	tridiagonal soma (t1.dim);
	
	for (int i=0; i<=t1.dim-2; i++) {
		soma.ld[i] = t1.ld[i]+t2.ld[i] ; soma.dd[i] = t1.dd[i]+t2.dd[i];
		soma.ud[i] = t1.ud[i]+t2.ud[i];
	}
	soma.dd[t1.dim-1] = t1.dd[t1.dim-1] + t2.dd[t1.dim-1];

	return soma;
}	

tridiagonal operator-(const tridiagonal &t1, const tridiagonal &t2)
{
	tridiagonal soma (t1.dim);
	
	for (int i=0; i<=t1.dim-2; i++) {
		soma.ld[i] = t1.ld[i]-t2.ld[i] ; soma.dd[i] = t1.dd[i]-t2.dd[i];
		soma.ud[i] = t1.ud[i]-t2.ud[i];
	}
	soma.dd[t1.dim-1] = t1.dd[t1.dim-1] - t2.dd[t1.dim-1];
	
	return soma;
}	

tridiagonal operator*(const Complex CC, const tridiagonal &tri)
{
	tridiagonal prod (tri.dim);
	for (int i=0; i<=tri.dim-2; i++) {
		prod.ld[i] = CC*tri.ld[i]; prod.dd[i] = CC*tri.dd[i];
		prod.ud[i] = CC*tri.ud[i];
	}
	prod.ld[tri.dim-1] = CC*tri.ld[tri.dim-1];
	return prod;
}

	
	
	

