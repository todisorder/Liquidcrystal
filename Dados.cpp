/*
 *  Dados.cpp
 *  
 *
 *  Created by Paulo Amorim on 12/Jul/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "Dados.h"

/********************************************************************/
//			Construtor vazio
/********************************************************************/

Dados::Dados ()
{
}

/********************************************************************/
//		Inicialização do parametro dadosdefinidos que evita 
//		a repetição da inserção de dados no caso de option = 1
//		em baixo.
/********************************************************************/

int Dados::dadosdefinidos=0;

/********************************************************************/
//		Função para definir os parâmetros do problema
/********************************************************************/

void Dados::set (void)
{

	/********************************************************************/
	/*		option = 0 --> parametros definidos manualmente.
	/*		option = 1 --> inserir numpt, dt, numiter.
	/********************************************************************/
	int option = 0;

	if (option == 0) {
		tfinal = 0.4;
		numpt =1000;
//		numiter = 200;
		dt = 0.00005;
		numiter = tfinal/dt;
		xmin = -15.;
		xmax = 15.;				//Calculado implicitamente
		dx = (xmax - xmin)/numpt;		
		delta = .05*dx;
		eps = .05*dx;		//viscosidade
        aa = 10.;
        HH = 10.;
        alpha = .3;
        beta = 0.29;
        gamma = .1;
        bb = 1.;
	}
	if (option == 1 && dadosdefinidos == 0) {
		cout << "numpt? "<< endl;
		cin >> numpt;
		cout << "dt? "<< endl;
		cin >> dt;
		cout << "numiter? "<< endl;
		cin >> numiter ;
		dadosdefinidos = 1;			//Para que não pergunte outra vez os dados.
	}
		focus = 1.;		//focus = 1 ou focus = -1; o problema certo é com = 1.		
//	dx = 0.1;
//	numpt = (xmax - xmin)/dx;
//		tfinal = numiter*dt;			//Calculado implicitamente
		BVLeftuu = 0.;		//Não uso BVs nesta versão
		BVRightuu = 0.;		//Não uso BVs nesta versão
		BVLeftvv = 0.;		//Não uso BVs nesta versão
		BVRightvv = 0.;		//Não uso BVs nesta versão

}

/********************************************************************/
//		Função para escrever os parametros 
/********************************************************************/

void Dados::write (void)
{
	cout << endl << "Dados:  \n";
	cout << " numpt = " << numpt << "\t"
    << " xmin = " << xmin << "\t" << " xmax = " << xmax << "\t" << " dx = " << dx << endl;
	cout << " numiter = " << numiter << "\t"
    << " tfinal = " << tfinal << "\t" << " dt = " << dt << endl << endl;
}

/********************************************************************/
//		Função para ir buscar os parametros int
/********************************************************************/

int Dados::get_int (const char *valor_inserido)
{
	if (!strncmp (valor_inserido, "numpt", 100))
    {
		return numpt;
    }
	else if (!strncmp (valor_inserido, "numiter", 100))
    {
		return numiter;
    }
	else
    {
		cout << "invalid request in get_int(). Exiting to system" << endl;
		exit (1);
    }
}

/********************************************************************/
//		Função para ir buscar os parametros double
/********************************************************************/

double Dados::get_double (const char *valor_inserido)
{
	if (!strncmp (valor_inserido, "dx", 100))
    {
		return dx;
    }
	else if (!strncmp (valor_inserido, "dt", 100))
    {
		return dt;
    }
	else if (!strncmp (valor_inserido, "xmin", 100))
    {
		return xmin;
    }
	else if (!strncmp (valor_inserido, "xmax", 100))
    {
		return xmax;
    }
	else if (!strncmp (valor_inserido, "tfinal", 100))
    {
		return tfinal;
    }
	else if (!strncmp (valor_inserido, "focus", 100))
    {
		return focus;
    }
	else if (!strncmp (valor_inserido, "eps", 100))
    {
		return eps;
    }	
    else if (!strncmp (valor_inserido, "delta", 100))
    {
        return delta;
    }
    else if (!strncmp (valor_inserido, "aa", 100))
    {
        return aa;
    }
    else if (!strncmp (valor_inserido, "HH", 100))
    {
        return HH;
    }
    else if (!strncmp (valor_inserido, "alpha", 100))
    {
        return alpha;
    }
    else if (!strncmp (valor_inserido, "beta", 100))
    {
        return beta;
    }
    else if (!strncmp (valor_inserido, "gamma", 100))
    {
        return gamma;
    }
    else if (!strncmp (valor_inserido, "bb", 100))
    {
        return bb;
    }
	else
    {
		cout << "invalid request in get_double(). Exiting to system" << endl;
		exit (1);
    }
}

/********************************************************************/
//		Função para ir buscar os parametros Complex
/********************************************************************/

Complex Dados::get_Complex (const char *valor_inserido)
{
	if (!strncmp (valor_inserido, "BVLeftuu", 100))
    {
		return BVLeftuu;
    }
	else if (!strncmp (valor_inserido, "BVRightuu", 100))
    {
		return BVRightuu;
    }
	else if (!strncmp (valor_inserido, "BVRightvv", 100))
    {
		return BVRightvv;
    }
	else if (!strncmp (valor_inserido, "BVLeftvv", 100))
    {
		return BVLeftvv;
    }
	else
    {
		cout << "invalid request in get_Complex(). Exiting to system" << endl;
		exit (1);
    }
}
