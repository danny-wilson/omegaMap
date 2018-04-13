objects = main.o omegaMap.o io.o mcmc.o tools.o
myutils = argumentwizard.h controlwizard.h DNA.h lotri_matrix.h \
	matrix.h MLST.h mutation.h mydouble.h myerror.h myutils.h \
	random.h utils.h vector.h
headers = omegaMap.h $(myutils)
idir = "./"
OLEVEL = 2

## INSTRUCT MAKE TO MAKE MULTIPLE GOALS
.PHONY : all
all : omegaMap omegaMapTP summarize decode order permute

## LINK
omegaMap : likelihood.o $(objects)
	g++ -w -O3 -o omegaMap likelihood.o $(objects)

omegaMapTP : likelihood.TP.o $(objects)
	-g++ -w -O3 -o omegaMapTP likelihood.TP.o $(objects)

summarize : summarize.o
	-g++ -w -O3 -o summarize summarize.o

decode : decode.o likelihood.o omegaMap.o io.o mcmc.o tools.o
	-g++ -w -O3 -o decode decode.o likelihood.o omegaMap.o io.o mcmc.o tools.o

order : order.o
	-g++ -w -O3 -o order order.o

permute: permute.o
	-g++ -w -O3 -o permute permute.o

## COMPILE
main.o : main.cpp $(headers)
	g++ -w -O3 -c -o main.o main.cpp -I$(idir)

omegaMap.o : omegaMap.cpp $(headers)
	g++ -w -O3 -c -o omegaMap.o omegaMap.cpp -I$(idir)

io.o : io.cpp $(headers)
	g++ -w -O3 -c -o io.o io.cpp -I$(idir)

likelihood.o : likelihood.cpp $(headers)
	g++ -w -O3 -c -o likelihood.o likelihood.cpp -I$(idir)

likelihood.TP.o : likelihood.cpp $(headers)
	-g++ -w -O3 -c -o likelihood.TP.o likelihood.cpp -I$(idir) -D_TESTPRIOR

mcmc.o : mcmc.cpp $(headers)
	g++ -w -O$(OLEVEL) -c -o mcmc.o mcmc.cpp -I$(idir)

tools.o : tools.c paml.h
	gcc -w -O3 -c -o tools.o tools.c -I$(idir)

summarize.o : summarize.cpp $(headers)
	-g++ -w -O3 -c -o summarize.o summarize.cpp -I$(idir)

decode.o : decode.cpp decode.h $(headers)
	-g++ -w -O3 -c -o decode.o decode.cpp -I$(idir)

order.o : order.cpp myerror.h random.h vector.h utils.h
	-g++ -w -O3 -c -o order.o order.cpp -I$(idir)

permute.o : permute.cpp $(myutils)
	-g++ -w -O3 -c -o permute.o permute.cpp -I$(idir)

## MAKE CLEAN
clean :
	-rm $(objects) likelihood.o likelihood.TP.o summarize.o decode.o order.o permute.o
