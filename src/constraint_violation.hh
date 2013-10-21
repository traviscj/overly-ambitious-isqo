
#ifndef __GUARD_CONSTRAINT_VIOLATION
#define __GUARD_CONSTRAINT_VIOLATION

#include <vector>

#include "step.hh"
#include "iterate.hh"
#include "nlp.hh"
#include "nlp_state.hh"


//! \brief A class returning a function for evaluating constraint violation of an iterate or step.
//!
//! Here, we evaluate the constraint violation 
//!   \f[ v(x):= ||c_{\mathcal{E}}(x)||_1 + ||[c_{\mathcal{I}}(x)]^+||_1 \f]
//! where \f$[x]^+ = max(x,0)\f$, or the constraint violation after a step,
//!   \f[ v(x,d):= ||c_{\mathcal{E}}(x) + J_{\mathcal{E}}\cdot d||_1 + ||[c_{\mathcal{I}}(x) + J_{\mathcal{I}}\cdot d]^+||_1 \f]
//! 
//! As always, we assume that the Nlp object returns a model of the form (for a specified \f$x\f$, \f$\lambda_{\mathcal{E}}\f$, and \f$\lambda_{\mathcal{I}}\f$)
//! \f{align*}{\min& \qquad  f(x) + \nabla f(x)^T d + \frac{1}{2}d^T\left[ \nabla^2_{xx} f(x) + \sum_{i\in\mathcal{E}} \nabla^2_{xx} c_{i}(x) \cdot \lambda_{(i)} + \sum_{i\in\mathcal{I}} \nabla^2_{xx} c_{i}(x) \cdot \lambda_{(i)} \right]d \\ \text{s.t.}& \qquad  J_{\mathcal{E}}\cdot d + c_{\mathcal{E}}(x) = 0\\& \qquad   J_{\mathcal{I}}\cdot d + c_{\mathcal{I}}(x) \leq 0 \f}
class ConstraintViolationFunction : public FunctionWithNLPState {
public:
    
    //! \brief construct a ConstraintViolationFunction
    //!
    //! \param nlp an Nlp object capable of returning constraint values and jacobians.
	ConstraintViolationFunction(iSQOControlPanel &control, Nlp &nlp);
    
    //! \brief evaluate the constraint violation function at a particular iterate
	double operator()(const iSQOIterate &iterate) const;
    //! \brief evaluate the constraint violation function using the linearized information from a particular iterate.
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) const;	
	
protected:
    //! \brief combine equality and inequality vectors from either of the operator() functions
    //! \param eq_items a vector of equality constraint values
    //! \param ieq_items a vector of inequality constraint values
    double two_vectors(const std::vector<double> &eq_items, const std::vector<double> &ieq_items) const;
};

#endif
