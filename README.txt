Instruções para 
/Users/Paulo/Matematica/Programas/liquidcrystal
2022
Paulo Amorim
/********************************************************************/
//		Resolve	o sistema de tipo Benney
//
//		iu_t +  u_xx = -ru + au|u|^2 + H^2 x^2 u	(NLS)
//		r_t = w
//		w_t - eps w_xx = (sig(r_x))_x - br + |u|^2	(liquid crystal)
//
//		com:
//		  sig(v) = bv + lambda v^3,
//		  lambda = (2/3)gamma*(alpha-beta).
//
//		com método semi-implícito Cranck-Nicholson + Newton para a
//		não-linearidade de NLS (?) e (?)
//
//		Paulo Amorim 2022.
/********************************************************************/

--> Para compilar, fazer 
>make visc
(ver o makefile)

--> Para executar, fazer 
>./exec

--> A execução cria ficheiros de dados com paths especificados em main.cpp.
(ver a parte "Nomes de Folders e files" em main.cpp)

--> Para visualizar com gnuplot:
>gnuplot
>cd "[path dos ficheiros]"
>plot "SolUURe.txt" w l, "SolUUAbs.txt" w l, "SolUUIm.txt" w l, "SolVV.txt" w l
(por exemplo). Ou então
>load plots.plt
ou outro ficheiro .plt.

ESTRUTURA dos ficheiros:

*************** main.cpp ******************
*************************************************
- Define instancias das classes Dados e Solver_SCL.
- Escreve dados iniciais
- Main loop:
		Prob.IteraClawLFImplicit(data);			//L-Friedrichs semi-implicito (OK)
		Prob.Newton(data);						//NLS acoplada com CLaw
- Escrever resultados em ficheiros.
*************************************************
*************************************************

*************** Solver_SCL.cpp ******************
*************************************************
--> construtor Solver_SCL::Solver_SCL(Dados *data)
- Inicializa tudo com uns Dados* data. Isto é muito importante (cf. linha 21).
- Define dados iniciais UU e VV - comentar o que não interessa
- Define Solução exacta UU e VV quando há - comentar o que não interessa
--> FIM DE construtor Solver_SCL::Solver_SCL(Dados *data)
--> Método de Newton para Nonlinear Schroedinger (NLS)
--> Solver NLS linear para testes
--> Def. da função fluxo f da Claw
--> Def. da função de acoplamento g e g' da Claw
--> Fluxos numéricos para Claw
--> 2 Esquemas explícitos INSTAVEIS
--> Esquema Lax-Friedrichs Implícito BOM mas não usa o fluxo numérico.
--> FIM DE Solver_SCL.cpp
*************************************************
*************************************************

*************** Dados.cpp ******************
*************************************************
--> Define os parametros do problema, tipo dx, dt, numpts, etc.
*************************************************
*************************************************

*************** Vectores.cpp ******************
*************************************************
Define a classe "vect".
Não esquecer:
-- Os destrutores são essenciais
-- Os operadores = e o copy constructor tem de ser feitos com ciclos tipo comps[i] = oldvector.comps[i],
	e não simplesmente comps = oldvector.comps.
-- Contém o produto de matriz TRIDIAGONAL com vector e resolução de sistema tridiagonal.
*************************************************
*************************************************

*************** Tridiagonal.cpp ******************
*************************************************
Define a classe "tridiagonal".
Atenção que as três diagonais são vectores normais, não são da classe vect !!
*************************************************
*************************************************