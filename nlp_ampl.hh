
#ifndef __GUARD_NLP_AMPL_HH
#define __GUARD_NLP_AMPL_HH

#include <vector>
#include <string>

#include <qpOASES.hpp>
// annoyingly... AMPL redefines a bunch of things, including 'filename', so it must come after the qpOASES header.
#include "asl.h"

#include "iterate.hh"
#include "nlp.hh"

class AmplNlp : public Nlp {
public:
	// AmplNlp(char *stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
	AmplNlp(std::string stub_str);
	
	~AmplNlp();
	
	iSQOIterate initial();
    
	// zeroth order NLP quantities:
	double objective(const iSQOIterate &iterate);
	
	std::vector<double> constraints_equality(const iSQOIterate &iterate);
	std::vector<double> constraints_inequality(const iSQOIterate &iterate);
	
	// first order NLP quantities:
	std::vector<double> objective_gradient(const iSQOIterate &iterate);
	matrix constraints_equality_jacobian(const iSQOIterate &iterate);
	qpOASES::SparseMatrix constraints_equality_jacobian_sparse(const iSQOIterate &iterate);
	matrix constraints_inequality_jacobian(const iSQOIterate &iterate);
	
	// second order NLP quantities:
	matrix lagrangian_hessian(const iSQOIterate &iterate);
	
	
	std::size_t num_lower_ieq();
	std::size_t num_upper_ieq();
    
private:
protected:
	bool PRINT_;
	ASL *asl_;
	FILE *nl_;
	
	std::vector<double> X0_;
	fint *nerror_;
	std::vector<size_t> equality_constraints_;
	std::vector<size_t> inequality_constraints_lower_;
	std::vector<size_t> inequality_constraints_upper_;
	
	std::vector<size_t> variable_equality_;
	std::vector<size_t> variable_bound_lower_;
	std::vector<size_t> variable_bound_upper_;
};

#endif
