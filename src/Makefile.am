UNITTEST_LIB=-L/usr/local/lib -lUnitTest++


lib_LIBRARIES = liboverly_ambitious_isqo-0.1.a
bin_PROGRAMS = test_matrices overly_ambitious_isqo

liboverly_ambitious_isqo_0_1_adir = .
liboverly_ambitious_isqo_0_1_a_SOURCES = constraint_violation.cc hessian_shifter.cc         isqo_config.cc iterate.cc line_search.cc linear_model.cc linear_model_reduction.cc matrix.cc nlp.cc nlp_ampl.cc nlp_hs014.cc              penalty_function.cc residual_function.cc solve_subproblem.cc step.cc subproblem.cc text_output.cc utilities.cc
liboverly_ambitious_isqo_0_1_a_HEADERS = constraint_violation.hh hessian_shifter.hh isqo.hh isqo_config.hh iterate.hh line_search.hh linear_model.hh linear_model_reduction.hh matrix.hh nlp.hh nlp_ampl.hh nlp_hs014.hh nlp_state.hh penalty_function.hh residual_function.hh solve_subproblem.hh step.hh subproblem.hh text_output.hh utilities.hh



overly_ambitious_isqodir = .
# TODO fix this:
overly_ambitious_isqo_LDADD=liboverly_ambitious_isqo-0.1.a /Users/traviscj/LocalPrograms/ampl/solvers/amplsolver.a
overly_ambitious_isqo_SOURCES = isqo_functor.cc 
overly_ambitious_isqo_HEADERS = isqo.hh

test_matricesdir = .
test_matrices_LDADD=liboverly_ambitious_isqo-0.1.a
test_matrices_SOURCES = test_matrices.cc 
test_matrices_HEADERS = matrix.hh