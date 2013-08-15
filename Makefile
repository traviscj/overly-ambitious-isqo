
# QPOASES_PATH=/Users/traviscj/optimization/qpOASES
QPOASES_PATH=/Users/traviscj/LocalPrograms/qpOASES-3.0beta
CFLAGS=-I${QPOASES_PATH}/include
LDFLAGS=-L${QPOASES_PATH}/bin -lqpOASES
# -lqpOASESextras
# ${QPOASES_PATH}/src/libqpOASES.a ${QPOASES_PATH}/src/libqpOASESextras.a
LDFLAGS+=-L/usr/local/Cellar/gfortran/4.8.1/gfortran/lib -lgfortran  -llapack -lblas

default: build run
	
run:
	./isqo_functor
build:
	g++ isqo_functor.cc ${CFLAGS} ${LDFLAGS} -o isqo_functor
	