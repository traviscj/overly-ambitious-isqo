
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
    
    void ConstructHelper(std::string stub_str);
	
	iSQOIterate initial();
    
	// zeroth order NLP quantities:
	double objective(const iSQOIterate &iterate);
	
	std::vector<double> constraints_equality(const iSQOIterate &iterate);
	std::vector<double> constraints_inequality(const iSQOIterate &iterate);
	
	// first order NLP quantities:
	std::vector<double> objective_gradient(const iSQOIterate &iterate);
	    
	virtual std::shared_ptr<matrix_base_class> constraints_equality_jacobian(const iSQOIterate &iterate) = 0;
	virtual std::shared_ptr<matrix_base_class> constraints_inequality_jacobian(const iSQOIterate &iterate) = 0;
	virtual std::shared_ptr<matrix_base_class> lagrangian_hessian(const iSQOIterate &iterate) = 0;
	
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
    void jacobian_update(const iSQOIterate &iterate);
    void hessian_update(const iSQOIterate &iterate);

    // TODO make these returned, let a caching class worry about it.
    // actually, not sure that this is possible--We have either of the constraints_*_jacobian_sparse functions call
    // sparse_jacobian_update, which in turn calls the jacobian evaluation functions.
    // the purpose of this is to NOT evaluate the constraints multiple times.
    // maybe we can only get away with not storing the hessian?
    // or maybe we should be more careful about which functions we 'coneval', which might help?
    // or... maybe we just throw caution to the wind for now, and clean it up later.
    std::shared_ptr<matrix_base_class> eq_jacobian_;
    std::shared_ptr<matrix_base_class> ieq_jacobian_;
    std::shared_ptr<matrix_base_class> hessian_;  
};

class DenseAmplNlp : public AmplNlp {
public:
	DenseAmplNlp(std::string stub_str);
	~DenseAmplNlp();

    std::shared_ptr<matrix_base_class> constraints_equality_jacobian(const iSQOIterate &iterate);
    std::shared_ptr<matrix_base_class> constraints_inequality_jacobian(const iSQOIterate &iterate);
    std::shared_ptr<matrix_base_class> lagrangian_hessian(const iSQOIterate &iterate);
private:
protected:

};

class SparseAmplNlp : public AmplNlp {
public:
	SparseAmplNlp(std::string stub_str);
	~SparseAmplNlp();
    
    std::shared_ptr<matrix_base_class> constraints_equality_jacobian(const iSQOIterate &iterate);
    std::shared_ptr<matrix_base_class> constraints_inequality_jacobian(const iSQOIterate &iterate);
    std::shared_ptr<matrix_base_class> lagrangian_hessian(const iSQOIterate &iterate);
private:
protected:
};

#endif
