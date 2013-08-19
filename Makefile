
# QPOASES_PATH=/Users/traviscj/optimization/qpOASES
QPOASES_PATH=/Users/traviscj/LocalPrograms/qpOASES-3.0beta

AMPL_LIB=/Users/traviscj/LocalPrograms/ampl/solvers/amplsolver.a
AMPL_INC=-I/Users/traviscj/LocalPrograms/ampl/solvers/

UNITTEST_LIB=-L/usr/local/lib -lUnitTest++

CFLAGS=-I${QPOASES_PATH}/include ${AMPL_INC}
LDFLAGS=-L${QPOASES_PATH}/bin -lqpOASES ${UNITTEST_LIB}
# -lqpOASESextras
# ${QPOASES_PATH}/src/libqpOASES.a ${QPOASES_PATH}/src/libqpOASESextras.a
LDFLAGS+=-L/usr/local/Cellar/gfortran/4.8.1/gfortran/lib -lgfortran  -llapack -lblas ${AMPL_LIB}

# test:
# 	g++ tester.cc ${CFLAGS} ${LDFLAGS} -o tester -O0 -g
# 	./tester

default: build run
	
run:
	./isqo_functor
build:
	g++ isqo_functor.cc ${CFLAGS} ${LDFLAGS} -o isqo_functor -O0 -g
	