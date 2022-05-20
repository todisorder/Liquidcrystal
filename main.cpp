/********************************************************************/
//		Resolve	a Equação visoelasticidade + Schrodinger
//
//
//
//
//		
//		COM RUNGE-KUTTA 4 EM TEMPO, METODO EXPLICITO
//		
//		Paulo Amorim 2010/2011.
/********************************************************************/


#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <string>
#include <algorithm> //Para max, min

#include "Dados.h"
#include "Tridiagonal.h"
#include "Vectores.h"
#include "Solver_Visc.h"


using namespace std;

const double pi = 3.14159265358979323846;



/********************************************************************/
//					MAIN
/********************************************************************/
int main () {

	/********************************************************************/
	//			Definir os dados do problema.
	/********************************************************************/
	Dados *data = new Dados;	
	data->set ();
	data->write ();
	int numpt = data->get_int("numpt");
	int numiter = data->get_int("numiter");
	double xmin = data->get_double("xmin");
	double xmax = data->get_double("xmax");
	double dx = data->get_double("dx");
	double dt = data->get_double("dt");
	double focus = data->get_double("focus");
	double eps = data->get_double("eps");	
	double NLinftyAcumulada = 0.;
	double MaxDifTemp = 0.;
	double tcurrent = 0.;
	vect qA(numpt);
	int Imax = 1;
	int ImaxL2 = 1;	
	double NL2Temp = 0.;
	double NL2 = 0.;
	
	/********************************************************************/
	//			Definir uma classe Solver_Visc
	//			Ver Solver_Visc.cpp
	/********************************************************************/
	Solver_Visc Prob(data);

	cout << "Norma L2 de UU = "<< NormaL2(Prob.UU,dx)<<endl;
	/********************************************************************/
	//		Escrever dados iniciais
	/********************************************************************/
    ofstream foutUUinire ("UU-Inire.txt");
    for (int j = 0; j <= numpt - 1; j++)
    {
        foutUUinire << real(Prob.XX.comps[j])<<"\t"<< real(Prob.UU.comps[j]) <<endl;
    }
    ofstream foutUUiniabs ("UU-Iniabs.txt");
    for (int j = 0; j <= numpt - 1; j++)
    {
        foutUUiniabs << real(Prob.XX.comps[j])<<"\t"<< abs(Prob.UU.comps[j]) <<endl;
    }
    ofstream foutRRini ("RR-Ini.txt");
    for (int j = 0; j <= numpt - 1; j++)
    {
        foutRRini << real(Prob.XX.comps[j])<<"\t"<< real(Prob.RR.comps[j]) <<endl;
    }
    ofstream foutWWini ("WW-Ini.txt");
    for (int j = 0; j <= numpt - 1; j++)
    {
        foutWWini << real(Prob.XX.comps[j])<<"\t"<< real(Prob.WW.comps[j]) <<endl;
    }
	/********************************************************************/
	//		Abrir ficheiro para grafico 3d
	/********************************************************************/
	ofstream fout3D1 ("Sol3DUURe.txt");
	ofstream fout3D11 ("Sol3DUUAbs.txt");
	ofstream fout3D2 ("Sol3DVV.txt");
	ofstream fout3D3 ("Sol3DWW.txt");
	
	/********************************************************************/
	//			Main loop
	/********************************************************************/
	for (int i=1; i<=numiter ; i++) {
		
		tcurrent = i*dt;
		
//        for (int j = 0; j <= numpt-1; j++)
//        {
//            Prob.NonlocalWW.comps[j] = Prob.NonlocalWW.comps[j] + dt*(exp(tcurrent-dt)*Prob.WW.comps[j]);
//        }
		////////////////////////////////|
		Prob.RKVISC(data, tcurrent);////|					//Iteration--------
		////////////////////////////////|
		
		if (i%50 == 0) {
			cout << "Iter restantes = "<<numiter - i<<endl;
		}
		
		//Actualizar o maximo do erro
		MaxDifTemp = MaxDif(Prob.UU, qA);
		if (MaxDifTemp > NLinftyAcumulada) {
			NLinftyAcumulada = MaxDifTemp;
			Imax = i;
		}
		//Acumulado da norma L2 do erro
		NL2Temp = NormaL2(Prob.UU - qA,dx);
		if (NL2Temp > NL2) {
			NL2 = NL2Temp;
			ImaxL2 = i;
		}
		
		/********************************************************************/
		//				Escrever a evolução
		/********************************************************************/
		int tresD = 1;		//0 --> no 3D; 1 --> 3D yes.		
		if (tresD == 1) {
			int samplt = numiter/min(numiter,40);		//Para ficar com x ptos em tempo
			int samplx = numpt/min(numpt,1000);			//Para ficar com x ptos em espaço
			if (i%samplt == 0) {	//Escrever no fich3d só de x em x iterações
                for (int j = 0; j <= numpt - 1; j=j+samplx)    //e só de y em y pontos espaciais
                {
                    fout3D1 <<  real(Prob.XX.comps[j])<<":";
                }
                fout3D1 << ";";
                for (int j = 0; j <= numpt - 1; j=j+samplx)    //e só de y em y pontos espaciais
                {
                    fout3D1 << real(Prob.UU.comps[j])<<":";
                }
                fout3D1 <<endl;
                
                for (int j = 0; j <= numpt - 1; j=j+samplx)    //e só de y em y pontos espaciais
                {
                    fout3D11 <<  real(Prob.XX.comps[j])<<":";
                }
                fout3D11 << ";";
                for (int j = 0; j <= numpt - 1; j=j+samplx)    //e só de y em y pontos espaciais
                {
                    fout3D11 << abs(Prob.UU.comps[j])<<":";
                }
                fout3D11 <<endl;
                
                for (int j = 0; j <= numpt - 1; j=j+samplx)    //e só de y em y pontos espaciais
                {
                    fout3D2 <<  real(Prob.XX.comps[j])<<":";
                }
                fout3D2 << ";";
                for (int j = 0; j <= numpt - 1; j=j+samplx)    //e só de y em y pontos espaciais
                {
                    fout3D2 << real(Prob.RR.comps[j])<<":";
                }
                fout3D2 <<endl;
                
                for (int j = 0; j <= numpt - 1; j=j+samplx)    //e só de y em y pontos espaciais
                {
                    fout3D3 <<  real(Prob.XX.comps[j])<<":";
                }
                fout3D3 << ";";
                for (int j = 0; j <= numpt - 1; j=j+samplx)    //e só de y em y pontos espaciais
                {
                    fout3D3 << real(Prob.WW.comps[j])<<":";
                }
                fout3D3 <<endl;
                
			}
		}

    
    
    
    
    }		//	Fim de iter
    
    
    
	NL2 = sqrt(NL2Temp);  //Norma L2 em (x,t)
	cout << "numiter = "<< numiter<<endl;	
	cout << "tfinal calculado = "<< dt*numiter<<endl;
	cout << "dt/dx = "<< dt/dx<<endl;
	cout << "dt/dx*dx = "<< dt/(dx*dx)<<endl;	
	cout << "Norma L2 de UU final = "<< NormaL2(Prob.UU,dx)<<endl;
	cout << "Focus = "<<focus<<endl;
	data->write ();
	/********************************************************************/
	//			Nomes de Folders e files
	//			ATENÇÃO: se as folders não existirem ele não as cria!!
	/********************************************************************/
	string FolderSols = "Sols";
	string FolderErrors = "Sols";	
	string stringSolUURe = FolderSols + "/SolUURe.txt";
	string stringSolUUIm = FolderSols + "/SolUUIm.txt";
	string stringSolUUAbs = FolderSols + "/SolUUAbs.txt";
	string stringSolUUExactRe = FolderSols + "/SolUUExactRe.txt";
	string stringSolUUExactIm = FolderSols + "/SolUUExactIm.txt";
	string stringSolUUExactAbs = FolderSols + "/SolUUExactAbs.txt";
	string stringErrUUL2 = FolderErrors + "/ErrUUL2.txt";
	string stringErrUULinfty = FolderErrors + "/ErrUULinfty.txt";
	string stringErrUUPontual = FolderErrors + "/ErrUUPontual.txt";	

	string stringSolVV = FolderSols + "/SolVV.txt";

	string stringSolWW = FolderSols + "/SolWW.txt";	
	
	char* FilenameSolUURe = new char[stringSolUURe.size()+1];
		strcpy (FilenameSolUURe,stringSolUURe.c_str());
	
//	const char* FilenameSolUUIm = (FolderSols + "/Sol-uu-Schr-Im.txt").c_str();
	char* FilenameSolUUIm = new char[stringSolUUIm.size()+1];
	strcpy (FilenameSolUUIm,stringSolUUIm.c_str());
	
//	const char* FilenameSolUUAbs = (FolderSols + "/Sol-uu-Schr-Abs.txt").c_str();
	char* FilenameSolUUAbs = new char[stringSolUUAbs.size()+1];
	strcpy (FilenameSolUUAbs,stringSolUUAbs.c_str());
	
//	const char* FilenameSolUUExactRe = (FolderSols + "/Sol-uu-Schr-EXACTAre.txt").c_str();
	char* FilenameSolUUExactRe = new char[stringSolUUExactRe.size()+1];
	strcpy (FilenameSolUUExactRe,stringSolUUExactRe.c_str());
	
//	const char* FilenameSolUUExactIm = (FolderSols + "/Sol-uu-Schr-EXACTAim.txt").c_str();
	char* FilenameSolUUExactIm = new char[stringSolUUExactIm.size()+1];
	strcpy (FilenameSolUUExactIm,stringSolUUExactIm.c_str());
	
//	const char* FilenameSolUUExactAbs = (FolderSols + "/Sol-uu-Schr-EXACTAabs.txt").c_str();
	char* FilenameSolUUExactAbs = new char[stringSolUUExactAbs.size()+1];
	strcpy (FilenameSolUUExactAbs,stringSolUUExactAbs.c_str());
		
//	const char* FilenameErrUUL2 = (FolderErrors + "/ErroL2UU.txt").c_str();
	char* FilenameErrUUL2 = new char[stringErrUUL2.size()+1];
	strcpy (FilenameErrUUL2,stringErrUUL2.c_str());

	char* FilenameErrUULinfty = new char[stringErrUULinfty.size()+1];
	strcpy (FilenameErrUULinfty,stringErrUULinfty.c_str());

	char* FilenameErrUUPontual = new char[stringErrUUPontual.size()+1];
	strcpy (FilenameErrUUPontual,stringErrUUPontual.c_str());
	
	char* FilenameSolVV = new char[stringSolVV.size()+1];
	strcpy (FilenameSolVV,stringSolVV.c_str());

	char* FilenameSolWW = new char[stringSolWW.size()+1];
	strcpy (FilenameSolWW,stringSolWW.c_str());

	
	/********************************************************************/
	//			Escrever resultados em ficheiros
	/********************************************************************/
	ofstream fout1 (FilenameSolUURe);
//    fout1 <<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
//    <<"# numpts = "<< numpt
//    << "\t Focus = " << focus << endl;
	for (int j = 0; j <= numpt - 1; j++)
    {
		fout1 << real(Prob.XX.comps[j])<<"\t"<< real(Prob.UU.comps[j]) <<endl;
    }
	ofstream fout2 (FilenameSolUUIm);
//    fout2 <<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
//    <<"# numpts = "<< numpt
//    << "\t Focus = " << focus << endl;
	for (int j = 0; j <= numpt - 1; j++)
    {
		fout2 << real(Prob.XX.comps[j])<<"\t"<< imag(Prob.UU.comps[j]) <<endl;
    }
	ofstream fout3 (FilenameSolUUAbs);
//    fout3 <<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
//    <<"# numpts = "<< numpt
//    << "\t Focus = " << focus << endl;
	for (int j = 0; j <= numpt - 1; j++)
    {
		fout3 << real(Prob.XX.comps[j])<<"\t"<< abs(Prob.UU.comps[j]) <<endl;
    }
	ofstream fout4 (FilenameSolVV);
//    fout4 <<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
//    <<"# numpts = "<< numpt
//    << "\t Focus = " << focus << endl;
	for (int j = 0; j <= numpt - 1; j++)
    {
		fout4 << real(Prob.XX.comps[j])<<"\t"<< real(Prob.RR.comps[j]) <<endl;
    }
	ofstream fout5 (FilenameSolWW);
//    fout5 <<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
//    <<"# numpts = "<< numpt
//    << "\t Focus = " << focus << endl;
	for (int j = 0; j <= numpt - 1; j++)
    {
		fout5 << real(Prob.XX.comps[j])<<"\t"<< real(Prob.WW.comps[j]) <<endl;
    }
	
	
	/********************************************************************/
	//			Escrever a Solução exacta (Quando há)
	/********************************************************************/
	ofstream foutUUexactaABS (FilenameSolUUExactAbs);	
	foutUUexactaABS <<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
	<<"# numpts = "<< numpt
	<< "\t Focus = " << focus << endl;
	for (int j = 0; j <= numpt - 1; j++)
	{
		foutUUexactaABS << real(Prob.XX.comps[j])<<"\t"<< abs(Prob.SOL_EXACTAUU.comps[j]) <<endl;
	}
	ofstream foutUUexactaRE (FilenameSolUUExactRe);	
	foutUUexactaRE <<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
	<<"# numpts = "<< numpt
	<< "\t Focus = " << focus << endl;
	for (int j = 0; j <= numpt - 1; j++)
	{
		foutUUexactaRE << real(Prob.XX.comps[j])<<"\t"<< real(Prob.SOL_EXACTAUU.comps[j]) <<endl;
	}
	ofstream foutUUexactaIM (FilenameSolUUExactIm);	
	foutUUexactaIM <<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
	<<"# numpts = "<< numpt
	<< "\t Focus = " << focus << endl;
	for (int j = 0; j <= numpt - 1; j++)
	{
		foutUUexactaIM << real(Prob.XX.comps[j])<<"\t"<< imag(Prob.SOL_EXACTAUU.comps[j]) <<endl;
	}
	ofstream foutErroPontual (FilenameErrUUPontual);	
	foutErroPontual <<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
	<<"# numpts = "<< numpt
	<< "\t Focus = " << focus << endl;
	for (int j = 0; j <= numpt - 1; j++)
	{
		foutErroPontual << real(Prob.XX.comps[j])<<"\t"<< abs(Prob.UU.comps[j] - Prob.SOL_EXACTAUU.comps[j]) <<endl;
	}
	
	/********************************************************************/
	//			Escrever algo em cada execução
	/********************************************************************/
	int history = 0; //
	if (history == 1) {
		ofstream ErroL2UU (FilenameErrUUL2, ios::app);
		ErroL2UU << "#=================================="<<endl
		<<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
		<<"# numpts = "<< numpt << "\t Iter do Max = " << ImaxL2 <<endl;
		ErroL2UU << "# epsilon, " <<"\t" << " Erro L^infty L^2 em (t,x)"<<endl;		
//		ErroL2UU << eps << "\t" << sqrt(NormaL2(Prob.UU - Prob.SOL_EXACTAUU, dx))<< endl;
		ErroL2UU << eps << "\t" << NL2 << endl;		
		//Erro L^\infty
		ofstream ErroLinfty (FilenameErrUULinfty, ios::app);
		ErroLinfty << "#=================================="<<endl
		<<"# dt = "<< dt<< "\t dx = "<<dx<< "\t Tfinal = "<<dt*numiter<<endl
		<<"# numpts = "<< numpt
		<< "\t xmin = " << xmin << "\t xmax = " << xmax << "\t Iter do Max = " << Imax<<endl;
		ErroLinfty << "# epsilon, " <<"\t" << " Erro L^infty em (t,x)"<<endl;
//		ErroLinfty << eps << "\t" << MaxDif(Prob.UU, Prob.SOL_EXACTAUU)<< endl;
		ErroLinfty << eps << "\t" << NLinftyAcumulada << endl;	
	}
	
	
    return 0;
}
/********************************************************************/
//			Fim de Main
/********************************************************************/
