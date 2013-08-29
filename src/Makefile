
QPOASES_PATH=/Users/traviscj/LocalPrograms/qpOASES-dj
QPOASES_INC=-I${QPOASES_PATH}/include_aw
QPOASES_LIB=-L${QPOASES_PATH}/bin -lqpOASES_aw

AMPL_PATH=/Users/traviscj/LocalPrograms/ampl/solvers
AMPL_INC=-I${AMPL_PATH}
AMPL_LIB=${AMPL_PATH}/amplsolver-debug.a

UNITTEST_LIB=-L/usr/local/lib -lUnitTest++

CFLAGS_INC=${QPOASES_INC} ${AMPL_INC} -I../include/
CFLAGS=-O0 -g -std=c++11 ${CFLAGS_INC}
LDFLAGS=${QPOASES_LIB} ${UNITTEST_LIB}
# -lqpOASESextras
# ${QPOASES_PATH}/src/libqpOASES.a ${QPOASES_PATH}/src/libqpOASESextras.a
LDFLAGS+=-L/usr/local/Cellar/gfortran/4.8.1/gfortran/lib -lgfortran  -llapack -lblas ${AMPL_LIB}

OBJECTS=utilities.o step.o iterate.o matrix.o nlp.o nlp_hs014.o nlp_ampl.o \
		constraint_violation.o penalty_function.o subproblem.o \
		residual_function.o solve_subproblem.o linear_model.o \
		linear_model_reduction.o line_search.o hessian_shifter.o \
		text_output.o

tests_src=$(wildcard test*.cc)
tests_bin=$(subst .cc,,${tests_src})

default: isqo-sparse isqo-dense

run: isqo
	./isqo

isqo-sparse: isqo_functor.cc ${OBJECTS}
	g++-4.8 -DISQO_USE_SPARSE ${OBJECTS} isqo_functor.cc -o isqo-sparse  ${CFLAGS} ${LDFLAGS}
isqo-dense: isqo_functor.cc ${OBJECTS}
	g++-4.8 ${OBJECTS} isqo_functor.cc -o isqo-dense  ${CFLAGS} ${LDFLAGS}

tests: ${tests_bin}

test_%: test_%.cc ${OBJECTS}
	g++-4.8 -o $@ $< ${OBJECTS} ${CFLAGS} ${LDFLAGS}

%.o: %.cc
	g++-4.8 -c $< ${CFLAGS}

valgrind: build
	valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./isqo_functor >> valgrind-output 2> valgrind-error

test:
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs001.nl > output/isqo_hs001.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs002.nl > output/isqo_hs002.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs003.nl > output/isqo_hs003.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs004.nl > output/isqo_hs004.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs005.nl > output/isqo_hs005.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs006.nl > output/isqo_hs006.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs007.nl > output/isqo_hs007.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs008.nl > output/isqo_hs008.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs009.nl > output/isqo_hs009.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs010.nl > output/isqo_hs010.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs011.nl > output/isqo_hs011.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs012.nl > output/isqo_hs012.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs013.nl > output/isqo_hs013.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs014.nl > output/isqo_hs014.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs015.nl > output/isqo_hs015.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs016.nl > output/isqo_hs016.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs017.nl > output/isqo_hs017.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs018.nl > output/isqo_hs018.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs019.nl > output/isqo_hs019.nl
	./isqo_functor ~/optimization/cute_nl_nopresolve/hs020.nl > output/isqo_hs020.nl

clean:
	@rm -f *.o
	@rm -f isqo ${tests_bin}
	@rm -rf *.dSYM
