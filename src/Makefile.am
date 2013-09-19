# QPOASES_PATH=/Users/traviscj/LocalPrograms/qpOASES-dj
# QPOASES_INC=-I${QPOASES_PATH}/include_aw
# QPOASES_LIB=-L${QPOASES_PATH}/bin -lqpOASES_aw

# AMPL_PATH=/Users/traviscj/LocalPrograms/ampl/solvers
# AMPL_INC=-I${AMPL_PATH}
# #AMPL_LIB=${AMPL_PATH}/amplsolver-debug.a
# #AMPL_LIB=${AMPL_PATH}/amplsolver-release.a
# AMPL_LIB=${AMPL_PATH}/amplsolver.a

UNITTEST_LIB=-L/usr/local/lib -lUnitTest++

# CFLAGS_INC=${QPOASES_INC} ${AMPL_INC}
# -I../include/
# AM_CPPFLAGS=-O3 -std=c++11 ${CFLAGS_INC}
# AM_LDFLAGS=${AM_LDFLAGS} ${AMPL_LIB}
#${QPOASES_LIB}
#${UNITTEST_LIB}
# -L/usr/local/Cellar/gfortran/4.8.1/gfortran/lib -lgfortran  -llapack -lblas 

overly_ambitious_isqo_includedir = .
overly_ambitious_isqodir = .
bin_PROGRAMS = overly_ambitious_isqo
overly_ambitious_isqo_LDADD=${ASL_LOCATION}/amplsolver.a
overly_ambitious_isqo_SOURCES = constraint_violation.cc hessian_shifter.cc isqo_config.cc isqo_functor.cc iterate.cc line_search.cc linear_model.cc linear_model_reduction.cc matrix.cc nlp.cc nlp_ampl.cc nlp_hs014.cc penalty_function.cc residual_function.cc solve_subproblem.cc step.cc subproblem.cc text_output.cc utilities.cc
overly_ambitious_isqo_HEADERS = constraint_violation.hh hessian_shifter.hh isqo.hh isqo_config.hh iterate.hh line_search.hh linear_model.hh linear_model_reduction.hh matrix.hh nlp.hh nlp_ampl.hh nlp_hs014.hh nlp_state.hh penalty_function.hh residual_function.hh solve_subproblem.hh step.hh subproblem.hh text_output.hh utilities.hh
