VISC_OBJS=main.o Dados.o  Vectores.o Solver_Visc.o
TESTES_OBJS=Testes.o Dados.o Vectores.o Solver_Visc.o 
CFLAGS=-Wall -g

visc:$(VISC_OBJS)
	g++ $(CFLAGS) -o exec $(VISC_OBJS) -lm

testes:$(TESTES_OBJS)
	g++ $(CFLAGS) -o testes $(TESTES_OBJS) -lm
clean:
	rm -f *.o;
	rm -f exec;

