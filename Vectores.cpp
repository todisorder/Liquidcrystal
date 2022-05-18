/*
 *  Vectores.cpp
 *  
 *
 *  Created by Paulo Amorim on 14/Jul/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include "Vectores.h"


//Constructors
//vect::vect()
//{
	
//	dim = 1;
//	comps = new Complex[dim];
//	for (int i=0; i<=dim-1; i++) {
//		comps[i] = 0.;
//	}
//}

vect::vect(Dados *data)		//Inicializar com Dados !!
{
//	data->set();		//não sei se é preciso...
	int nn = data->get_int("numpt");
	dim = nn;
	comps = new Complex[nn];
	for (int i=0; i<=dim-1; i++) {
		comps[i] = 0.;
	}	
}


vect::vect(int nn)
{
	dim = nn;
	comps = new Complex[nn];
	for (int i=0; i<=dim-1; i++) {
		comps[i] = 0.;
	}
}

vect::vect(int nn, Complex C)
{
	dim = nn;
	comps = new Complex[nn];
	for (int i=0; i<=dim-1; i++) {
		comps[i] = C ;
	}
}


vect::~vect()
{
	delete[] comps;
//	cout << "Hi"<<endl;			//Aparece!!!
}






//===========================================
//Copy constructor
vect::vect(const vect& u)
{
	dim = u.dim;
	comps = new Complex[dim];
	for (int i=0; i<=dim-1; i++) {
		comps[i] = u.comps[i] ;
	}
}

//===========================================
void vect::show() const	{
	for (int i=0; i<dim; i++) {
		cout << comps[i]<<endl;
	}
}
	
	

vect &vect::operator=(const vect &ve)
{
	//Vamos tentar isto::
	if (dim < ve.dim) {
		comps = new Complex[ve.dim];
	}
	dim=ve.dim;
	for(int i=0;i<dim;i++)
		comps[i]=ve.comps[i];
//	comps = ve.comps; //ASSIM NAO DA!!!! CUIDADO, MUDA OS DOIS!!!
	return *this;
}

vect &vect::operator+=(const vect &v2)
{
	dim=v2.dim;
    for(int i=0;i<dim;i++)
		comps[i]+=v2.comps[i];
	return *this;
}

vect &vect::operator-=(const vect &v2)
{
	dim=v2.dim;
    for(int i=0;i<dim;i++)
		comps[i]-=v2.comps[i];
	return *this;
}

vect vect::operator-( ) const
{
	vect vec(dim);
	for(int i=0;i<dim;i++)
		vec.comps[i]=-comps[i];
	return vec;
}

vect operator+(const vect &v1,const vect &v2)
{
	vect vec(v1.dim);
//	assert(v1.dim==v2.dim);
	for(int i=0;i<v1.dim;i++)
		vec.comps[i]=v1.comps[i]+v2.comps[i];
	return vec;
}

vect operator+(const vect &v1,const Complex ivalue )//vector plus a constant. 
{
	vect vec(v1.dim);
	
	for(int i=0;i<v1.dim;i++)
		vec.comps[i]=v1.comps[i]+ivalue;
	return vec;
}

vect operator-(const vect &v1,const Complex ivalue )//vector minus a constant. 
{
	vect vec(v1.dim);
	
	for(int i=0;i<v1.dim;i++)
		vec.comps[i]=v1.comps[i]-ivalue;
	return vec;
}


vect operator-(const vect &v1,const vect &v2)
{
	vect vec(v1.dim);
//	assert(v1.dim==v2.dim);
	for(int i=0;i<v1.dim;i++)
		vec.comps[i]=v1.comps[i]-v2.comps[i];
	return vec;
}


Complex operator*(const vect &v1,const vect &v2) //inner product between vectors
{
	Complex retrn=Complex(0.,0.);
//	assert(v1.dim==v2.dim);
    for(int i=0;i<v1.dim;i++)
		retrn+=v1.comps[i]*v2.comps[i];
	
	return retrn;
}


vect operator*(const Complex rho, const vect &vec)//constant vector product 
{
	vect ve(vec.dim);
	for(int i=0;i<vec.dim;i++)
		ve.comps[i]=rho*vec.comps[i]; 
	return ve;
}

vect operator*(const vect &vec,const tridiagonal &tri) //Tridiagonal matrix vector product
{
	vect result(vec.dim);
	result.comps[0] = tri.dd[0]*vec.comps[0] + tri.ud[0]*vec.comps[1];
	for (int i=1; i<=vec.dim-2; i++) {
		result.comps[i] = tri.ld[i-1]*vec.comps[i-1] + tri.dd[i]*vec.comps[i] + tri.ud[i]*vec.comps[i+1];
	}
	result.comps[vec.dim-1] = tri.ld[vec.dim-2]*vec.comps[vec.dim-2] + tri.dd[vec.dim-1]*vec.comps[vec.dim-1];	
	return result;
}

vect vect::realvect (void)
{
	vect result(dim);
	for(int i=0;i<dim;i++)
		result.comps[i]=real(comps[i]); 
	return result;
}

vect vect::conjugate (void)
{
	vect result(dim);
	for(int i=0;i<dim;i++)
		result.comps[i]=conj(comps[i]); 
	return result;
}

double MaxDif(const vect &u1, const vect &u2)
{
	int numpt = u1.dim;
	double max = abs(u1.comps[0] - u2.comps[0]);
	double max1 = 0.;
	for (int j=0; j<=numpt-1; j++) {
		max1 = abs(u1.comps[j] - u2.comps[j]);
		if (max1 > max ) {
			max = max1;
		}
	}
	return max;
}


vect vect::SolveTri(const tridiagonal& mat)
{		//Solves mat * x = uu
//	cout << "before ="<<endl;
//	(*this).show();
	int j;
	int nj = dim - 1;
	Complex beta;
	vect gamma(dim);
//	vect rhs(*this);
//	vect rhs(dim);
//	rhs.comps = (*this).comps;
	//aviso 1 aqui-----
	comps[0] = comps[0] / (beta = mat.dd[0]);
	for (j=1; j<=nj; j++) {
//		cout << "j="<<j<<endl;
		gamma.comps[j] = mat.ud[j-1] / beta;
		beta = mat.dd[j] - mat.ld[j-1] * gamma.comps[j];
		//Aviso 2 aqui------
		comps[j] = (comps[j] - mat.ld[j-1] * comps[j-1]) / beta;
	}
	for (j=(nj-1); j>=0; j--) {
//		cout << "jj="<<j<<endl;		
		comps[j] -= gamma.comps[j+1] * comps[j+1];
	}
//	cout<<"after="<<endl;
//	(*this).show();

	return *this;
	
//	mat.show();
//	cout<<"add this="<<comps[2]<<endl;
//	delete rhs.comps;
	
}

/*		Copia de SEGURANÇA
vect vect::SolveTri(const tridiagonal& mat)
{		//Solves mat * x = uu
	int j;
	int nj = dim - 1;
	Complex beta;
	vect gamma(dim);
	vect rhs(*this);
	//	vect rhs(dim);
	//	rhs.comps = (*this).comps;
	//aviso 1 aqui-----
	comps[0] = rhs.comps[0] / (beta = mat.dd[0]);
	for (j=1; j<=nj; j++) {
		//		cout << "j="<<j<<endl;
		gamma.comps[j] = mat.ud[j-1] / beta;
		beta = mat.dd[j] - mat.ld[j-1] * gamma.comps[j];
		//Aviso 2 aqui------
		comps[j] = (rhs.comps[j] - mat.ld[j-1] * comps[j-1]) / beta;
	}
	for (j=(nj-1); j>=0; j--) {
		//		cout << "jj="<<j<<endl;		
		comps[j] -= gamma.comps[j+1] * comps[j+1];
	}
	
	//	mat.show();
	rhs.show();cout<<"add="<<rhs.comps<<endl;cout<<"add this="<<this<<endl;
	//	delete rhs.comps;
	
}
*/

double NormaL2 (const vect& uu, double dx)
{
	double result = 0.;
	for (int i=0; i<uu.dim; i++) {
		result += dx*abs(uu.comps[i]*uu.comps[i]);
	}
	return result;
}

vect Compare  (const vect& u1, const vect& u2)
{
	return u1 - u2;
}




