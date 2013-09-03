
#ifndef __GUARD_NLP_AMPL_HH
#define __GUARD_NLP_AMPL_HH

#include <vector>
#include <string>

#include <qpOASES.hpp>
// annoyingly... AMPL redefines a bunch of things, including 'filename', so it must come after the qpOASES header.
#include "asl.h"

#include "iterate.hh"
#include "nlp.hh"
#include "utilities.hh"

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

class DenseAmplNlp : public AmplNlp {
public:
	// AmplNlp(char *stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
	DenseAmplNlp(std::string stub_str);
	
	~DenseAmplNlp() { }
	    
	matrix constraints_equality_jacobian(const iSQOIterate &iterate);
	matrix constraints_inequality_jacobian(const iSQOIterate &iterate);
	matrix lagrangian_hessian(const iSQOIterate &iterate);
	
	std::size_t num_lower_ieq();
	std::size_t num_upper_ieq();
    
    
	sparse_matrix constraints_equality_jacobian_sparse(const iSQOIterate &iterate) { throw "not implemented"; }
	sparse_matrix constraints_inequality_jacobian_sparse(const iSQOIterate &iterate) { throw "not implemented"; }
	sparse_matrix lagrangian_hessian_sparse(const iSQOIterate &iterate) { throw "not implemented"; }
    
private:
protected:
};

class SparseAmplNlp : public DenseAmplNlp {
public:
	// AmplNlp(char *stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
	SparseAmplNlp(std::string stub_str);
	
	~SparseAmplNlp() { }
	
	// first order NLP quantities:
	sparse_matrix constraints_equality_jacobian_sparse(const iSQOIterate &iterate);
	sparse_matrix constraints_inequality_jacobian_sparse(const iSQOIterate &iterate);
	// second order NLP quantities:
	sparse_matrix lagrangian_hessian_sparse(const iSQOIterate &iterate);
	
    
    // matrix constraints_equality_jacobian(const iSQOIterate &iterate) { throw "not implemented"; } 
    // matrix constraints_inequality_jacobian(const iSQOIterate &iterate) { throw "not implemented"; }
    // // second order NLP quantities:
    // matrix lagrangian_hessian(const iSQOIterate &iterate) { throw "not implemented"; }
	
private:
protected:
    void sparse_jacobian_update(const iSQOIterate &iterate);
    void sparse_hessian_update(const iSQOIterate &iterate);

    // TODO make these returned, let a caching class worry about it.
    // actually, not sure that this is possible--We have either of the constraints_*_jacobian_sparse functions call
    // sparse_jacobian_update, which in turn calls the jacobian evaluation functions.
    // the purpose of this is to NOT evaluate the constraints multiple times.
    // maybe we can only get away with not storing the hessian?
    // or maybe we should be more careful about which functions we 'coneval', which might help?
    // or... maybe we just throw caution to the wind for now, and clean it up later.
    sparse_matrix sparse_eq_jacobian_;
    sparse_matrix sparse_ieq_jacobian_;
    sparse_matrix sparse_hessian_;
};

#endif
