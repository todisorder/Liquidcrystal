/*
 *  Dados.h
 *  
 *
 *  Created by Paulo Amorim on 12/Jul/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _DADOS_H_
#define _DADOS_H_

#include <cmath>
#include <complex>

#include <iostream>
#include <string>

using namespace std;

typedef complex<double> Complex;

class   Dados
{
	
private:
	
	int     numpt, numiter;
	double  dt, dx, xmin, xmax, tfinal, focus, eps, Tzero, delta, aa, HH, alpha, beta, gamma, bb;
	Complex BVLeftuu;
	Complex BVLeftvv;
	Complex BVRightuu;
	Complex BVRightvv;
	static int dadosdefinidos;
	
public:
	
	Dados ();
	void    set (void);
	void    write (void);
	int     get_int (const char *valor_inserido);
	double  get_double (const char *valor_inserido);
	Complex get_Complex (const char *valor_inserido);
	
};

#endif /* ifndef _DADOS_H_ */
