
# QPOASES_PATH=/Users/traviscj/optimization/qpOASES
QPOASES_PATH=/Users/traviscj/LocalPrograms/qpOASES-3.0beta
AMPL_LIB=/Users/traviscj/LocalPrograms/ampl/solvers/amplsolver.a
AMPL_INC=-I/Users/traviscj/LocalPrograms/ampl/solvers/
CFLAGS=-I${QPOASES_PATH}/include ${AMPL_INC}
LDFLAGS=-L${QPOASES_PATH}/bin -lqpOASES
# -lqpOASESextras
# ${QPOASES_PATH}/src/libqpOASES.a ${QPOASES_PATH}/src/libqpOASESextras.a
LDFLAGS+=-L/usr/local/Cellar/gfortran/4.8.1/gfortran/lib -lgfortran  -llapack -lblas ${AMPL_LIB}

default: build run
	
run:
	./isqo_functor
build:
	g++ isqo_functor.cc ${CFLAGS} ${LDFLAGS} -o isqo_functor
	