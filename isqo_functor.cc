// overly-ambitious-iSQO
//  - Travis C. Johnson (traviscj@traviscj.com)
//  - August 2013

// basic philosophy is:
// - keep main loop as clean and as close to algo in the paper as humanly possible, at the expense of anything else.
// - don't try to fit functions (e.g. linear reduction) into classes for iterates
//   - INSTEAD: write functions that stand on their own (but functors, so they can have some persistent state)
// - don't try to fit multiple types of iterates into one class
//   - INSTEAD: separate the iterates from steps from subproblems from hessian shifting from...
// - don't try to fit more functionality than absolutely necessary into classes.
//   - INSTEAD: have easy-to-debug classes with a single simple point: handling hessian shifts, for instance.
// - don't try to pass in the kitchen sink simply to be able to evaluate some basic things.
//   - INSTEAD: Pass the NLP class, which returns f(x), c(x), \bar{c}(x), take it from there.
// - DO NOT use short (1-2 character) or greek symbols as variable names.
//   - INSTEAD: spell out exactly what you are iterating over. Primal indices like "i"? Try primal_index. Linesearch step size called "alpha"? Try linesearch_step_size.
//   - Exception: External class interfaces are a bit messier.
//   - Gotcha: A bit more consistent naming scheme would improve things even more.
// - whenever possible, return ONE object by value. 
//   - INSTEAD: let C++11x worry about RVO. 
//   - Exception: iSQOStep class... but maybe it doesn't need to be?
// - whenever possible, pass required objects by const reference.
// - whenever possible, use std::vector<double> instead of any other std::vectors.
// - for unfinished features, use assert to ensure they are not used.

// future feature list:
// - DONE !!! variable bounds !!!
// - DONE !!! more robust shifting strategy and/or Quasi-Newton. (Quasi-Newton might still be an interesting experiment)
// - DONE... rewrite some of the matrix computations into the matrix class. (But still need to check for stragglers...)
// - constraint/objective scaling (should fix HS99... )
// - other iSQO features (algo II-features.)
// - iQP interface! -- but how to do fallback/etc?
// - integrate DJ/CK/AW/et al's work on sparse qpOASES.
// - logbook-style reporting/logging, with sql.
// - use qpoases matrices instead of custom matrix code
// - sparse matrices from ampl/to qpOASES.
// - BQPD interface?
// - read parameters from files, store in some structure.
// - second order correction (might fix HS067...)
// - use multiple qp objects for getting subproblem solutions.
// - const all applicable functions.
// - 

// justification for weirder choices:
// - I use "functors" for some of the paper-defined functions, to capture program state. They work by overloading operator(), which is a big weird, but they are essentially functions with state. One key thing: State is still minimized in them--most are just wrappers around basic NLP functions. With this I was trying to accomplish clear input-output relations, instead of inventing functions that do incremental changes to a huge block of state. One downside to this, as written, is that many things get re-evaluated(This could be minimized by some caching of )

// future ideas:
// - would like to extend the FunctionWithNLPState class to cache the contents of the arguments, and return the last value if they haven't changed. I think this could be somewhat straightforward, but maybe I'm underestimating. This is another major "point" to doing 
// 

// helpful hints:
// - OSX: for alloc bugs, run with: export DYLD_INSERT_LIBRARIES=/usr/lib/libgmalloc.dylib. On other systems, valgrind.
// - currently, here, the reported KKT is the two-norm... in the matlab implementation, it is the inf norm, scaled by a bunch of things.

// remaining potential bugs:
// hs067.out:Failure - Did not converge - tiny step length... maybe bad step from QP solver? -- might be best fixed by a second-order correction.
// hs099.out:Failure - Did not converge - this *should* be fixed with problem rescaling - should add problem rescaling.
// other HS failures:
// hs087.out:Failure - Did not converge - was also iter limit in MATLAB implmentation.
// hs098.out:Failure - Did not converge - subproblem failure in MATLAB implementation.
// hs106.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs109.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs114.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs117.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs99exp.out:Failure - Did not converge - iteration limit in MATLAB implementation.

// standard libraries:
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>

// 
#include <qpOASES.hpp>

// iSQO headers:
#include "utilities.hh"
#include "step.hh"
#include "iterate.hh"
#include "nlp.hh"
#include "nlp_hs014.hh"
#include "nlp_ampl.hh"
#include "nlp_state.hh"
#include "constraint_violation.hh"

// using namespace std;std::std::vector

class PenaltyFunction : public FunctionWithNLPState{
public:
	PenaltyFunction(Nlp &nlp) : 
			FunctionWithNLPState(nlp), 
			constraint_violation_func_(nlp)
			{
		// std::cout << "-- Initializing penalty function." << std::endl;
	}
	double operator()(const iSQOIterate &iterate) const {
		// std::cout << "--- Calling the PenaltyFunction Functor..." << std::endl;
		double f = nlp_->objective(iterate);
		// std::cout << "--- objective: " << f << std::endl;
		return iterate.penalty_parameter_*f + constraint_violation_func_(iterate);
	}
protected:
	ConstraintViolationFunction constraint_violation_func_;
	// double penalty_parameter_;
};

class iSQOQuadraticSubproblem : public FunctionWithNLPState{
public:
	iSQOQuadraticSubproblem(Nlp &nlp, const iSQOIterate &iterate) :
				FunctionWithNLPState(nlp),
				num_qp_variables_(nlp.num_primal() + 2*nlp.num_dual()), 
				num_qp_constraints_(nlp.num_dual()),
				num_nlp_variables_(nlp.num_primal()), 
				num_nlp_constraints_eq_(nlp.num_dual_eq()), 
				num_nlp_constraints_ieq_(nlp.num_dual_ieq()),
				hessian_(num_qp_variables_, num_qp_variables_), 
				unshifted_hessian_diagonal_(num_qp_variables_), 
				first_shift_(false),
				jacobian_(num_qp_constraints_, num_qp_variables_), 
				gradient_(num_qp_variables_),
				lower_bound_(num_qp_variables_), 
				upper_bound_(num_qp_variables_), 
				jacobian_lower_bound_(num_qp_constraints_), 
				jacobian_upper_bound_(num_qp_constraints_),
				nlp_hessian_(num_nlp_variables_, num_nlp_variables_),
				nlp_eq_jacobian_(num_nlp_constraints_eq_, num_nlp_variables_),
				nlp_ieq_jacobian_(num_nlp_constraints_ieq_, num_nlp_variables_),
				nlp_objective_gradient_(num_nlp_variables_)
	 {
		nlp_objective_gradient_ = nlp_->objective_gradient(iterate);
		nlp_hessian_ = nlp_->lagrangian_hessian(iterate);
		nlp_eq_jacobian_ = nlp_->constraints_equality_jacobian(iterate);
		nlp_ieq_jacobian_ = nlp_->constraints_inequality_jacobian(iterate);
	
		for (size_t eq_constraint_index=0; eq_constraint_index < iterate.num_dual_eq_; ++eq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				jacobian_.set(eq_constraint_index,variables, nlp_eq_jacobian_.get(eq_constraint_index,variables));
			}
			jacobian_.set(eq_constraint_index,nlp_->num_primal()+eq_constraint_index, -1.0);
			jacobian_.set(eq_constraint_index,nlp_->num_primal()+iterate.num_dual_eq_+eq_constraint_index, +1.0);
		}
				
		for (size_t ieq_constraint_index=0; ieq_constraint_index < iterate.num_dual_ieq_; ++ieq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index,variables, nlp_ieq_jacobian_.get(ieq_constraint_index,variables));
			}
			jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index, nlp_->num_primal()+2*iterate.num_dual_eq_+ieq_constraint_index, -1.0);
			jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index, nlp_->num_primal()+2*iterate.num_dual_eq_+iterate.num_dual_ieq_+ieq_constraint_index, +1.0);
		}		
		for (size_t r=0; r<iterate.num_primal_; ++r) {
			for (size_t c=0; c<iterate.num_primal_; ++c) {
				hessian_.set(r,c, nlp_hessian_.get(r,c));
			}
		}
				
		// NLP gradient copied over
		for (size_t primal_index=0; primal_index<iterate.num_primal_; ++primal_index)
			gradient_[primal_index] = iterate.penalty_parameter_*nlp_objective_gradient_[primal_index];
		// equality slacks both have positive:
		for (size_t dual_eq_index=0; dual_eq_index<iterate.num_dual_eq_; ++dual_eq_index) {
			gradient_[iterate.num_primal_+dual_eq_index] = 1.0;
			gradient_[iterate.num_primal_+iterate.num_dual_eq_+dual_eq_index] = 1.0;
		}
		// inequality slacks only penalize the positive parts:
		for (size_t dual_ieq_index=0; dual_ieq_index<iterate.num_dual_ieq_; ++dual_ieq_index) {
			gradient_[iterate.num_primal_+2*iterate.num_dual_eq_+dual_ieq_index] = 1.0;
			gradient_[iterate.num_primal_+2*iterate.num_dual_eq_+iterate.num_dual_ieq_+dual_ieq_index] = 0.0;
		}
		
		std::vector<double> con_values_eq=nlp_->constraints_equality(iterate);
		std::vector<double> con_values_ieq=nlp_->constraints_inequality(iterate);
		
		for (size_t eq_index=0; eq_index<iterate.num_dual_eq_; ++eq_index) {
			jacobian_lower_bound_[eq_index] = -con_values_eq[eq_index];
			jacobian_upper_bound_[eq_index] = -con_values_eq[eq_index];
		}
		for (size_t ieq_index=0; ieq_index<iterate.num_dual_ieq_; ++ieq_index) {
			jacobian_lower_bound_[iterate.num_dual_eq_+ieq_index] = -INFINITY;
			jacobian_upper_bound_[iterate.num_dual_eq_+ieq_index] = -con_values_ieq[ieq_index];
		}
		for (size_t variable_index=0; variable_index < iterate.num_primal_; ++variable_index) {
			lower_bound_[variable_index] = -1e10;
			upper_bound_[variable_index] = +1e10;
		}
		for (size_t variable_index=0; variable_index < 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_); ++variable_index) {
			lower_bound_[iterate.num_primal_+variable_index] = 0.0;
			upper_bound_[iterate.num_primal_+variable_index] = 1e10;
		}
	}
	
	void inc_regularization(double hessian_shift) {
		bool PRINT = false;
		// make a copy of the diagonal entries, if we haven't already got one...
		if (! first_shift_) {
			if (PRINT) std::cout << "making the first-time-only copy..." << std::endl;
			for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
					unshifted_hessian_diagonal_[diagonal_entry] = hessian_.get(diagonal_entry,diagonal_entry);
			}
			first_shift_ = true;
		}
		
		for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
				hessian_.set(diagonal_entry,diagonal_entry, unshifted_hessian_diagonal_[diagonal_entry] + hessian_shift);
		}
		for (size_t diagonal_entry=num_nlp_variables_; diagonal_entry < num_qp_variables_; ++diagonal_entry) {
			hessian_.set(diagonal_entry, diagonal_entry, hessian_shift);
		}
	}
	
	void print() const {
		std::cout << std::endl;
		std::cout << "H = [";
		std::cout << hessian_;
		std::cout << "];" << std::endl;
		std::cout << "A = [";
		std::cout << jacobian_;
		std::cout << "];" << std::endl;
		
		std::cout << "g = " << gradient_ << "';" << std::endl;
		std::cout << "lb = " << lower_bound_ << "';" << std::endl;
		std::cout << "ub = " << upper_bound_ << "';" << std::endl;
		std::cout << "lbA = " << jacobian_lower_bound_ << "';" << std::endl;
		std::cout << "ubA = " << jacobian_upper_bound_ << "';" << std::endl;
	}
	int num_qp_variables_, num_qp_constraints_;
	int num_nlp_variables_, num_nlp_constraints_eq_, num_nlp_constraints_ieq_;
	matrix hessian_;
	bool first_shift_;
	std::vector<double> unshifted_hessian_diagonal_;
	matrix jacobian_;
	std::vector<double> gradient_;
	std::vector<double> lower_bound_;
	std::vector<double> upper_bound_;
	std::vector<double> jacobian_lower_bound_;
	std::vector<double> jacobian_upper_bound_;
	matrix nlp_hessian_;
	matrix nlp_eq_jacobian_;
	matrix nlp_ieq_jacobian_;
	std::vector<double> nlp_objective_gradient_;
private:
protected:
};

class ResidualFunction : public FunctionWithNLPState {
public:
	ResidualFunction(Nlp &nlp) : FunctionWithNLPState(nlp) {
		// std::cout << "-- Initializing nlp residual function." << std::endl;
	}
	double resid_helper(const iSQOIterate &iterate, std::vector<double> stationarity, std::vector<double> constraint_eq_values, std::vector<double> constraint_ieq_values, std::vector<double> constraint_eq_dual_values, std::vector<double> constraint_ieq_dual_values) const {
		std::vector<double> rho(stationarity.size() + 2*constraint_eq_values.size() + 2*constraint_ieq_values.size());
		
		matrix jac_ce = nlp_->constraints_equality_jacobian(iterate);
		matrix jac_ci = nlp_->constraints_inequality_jacobian(iterate);
				
		size_t rho_index = 0;

		double eq_signflip = 1;
		double ieq_signflip = 1;
		// std::cout << "sta pre: " << stationarity << std::endl << std::endl << std::endl;
		std::vector<double> eq_jacobian_trans_times_mults = jac_ce.multiply_transpose(constraint_eq_dual_values);
		std::vector<double> ieq_jacobian_trans_times_mults = jac_ci.multiply_transpose(constraint_ieq_dual_values);
		
		for (size_t stationarity_index=0; stationarity_index < stationarity.size(); ++stationarity_index) {
			rho[rho_index] = stationarity[stationarity_index] + eq_jacobian_trans_times_mults[stationarity_index] + ieq_jacobian_trans_times_mults[stationarity_index];
			++rho_index;
		}
		for (size_t constraint_eq_index=0; constraint_eq_index < nlp_->num_dual_eq(); ++constraint_eq_index) {
			rho[rho_index] = std::min(bracket_plus(constraint_eq_values[constraint_eq_index]), 1.0 - eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
			++rho_index;
			
			rho[rho_index] = std::min(bracket_minus(constraint_eq_values[constraint_eq_index]), 1.0 + eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
			++rho_index;
		}
		
		for (size_t constraint_ieq_index=0; constraint_ieq_index < nlp_->num_dual_ieq(); ++constraint_ieq_index) {
			rho[rho_index] = std::min(bracket_plus(constraint_ieq_values[constraint_ieq_index]), 1.0 - ieq_signflip*constraint_ieq_dual_values[constraint_ieq_index]);
			++rho_index;

			rho[rho_index] = std::min(bracket_minus(constraint_ieq_values[constraint_ieq_index]), 0.0 + ieq_signflip*constraint_ieq_dual_values[constraint_ieq_index]);
			++rho_index;
			
		}
		// std::cout << "rho(" << iterate.penalty_parameter_ << "): " << rho << std::endl;
		double total = 0.0;
		for (size_t rho_index = 0; rho_index < rho.size(); ++rho_index)
			total += rho[rho_index]*rho[rho_index];
		
		return sqrt(total);
	}
	double operator()(const iSQOIterate &iterate) const {
		bool PRINT=false;
		// std::cout << "OPERATOR FOR RESID FUNC: " << iterate.penalty_parameter_ << std::endl;
		
		std::vector<double> stationarity = nlp_->objective_gradient(iterate);
		
		for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
			// mu*grad f:
			stationarity[stationarity_index] *= iterate.penalty_parameter_;
		}
		
		std::vector<double> eq_con_values = nlp_->constraints_equality(iterate);
		std::vector<double> ieq_con_values = nlp_->constraints_inequality(iterate);
		
		if (PRINT) {
			std::cout << std::endl << "mu: " << iterate.penalty_parameter_ << std::endl;
			std::cout << "sta (pre): " << nlp_->objective_gradient(iterate) << std::endl;
			std::cout << "sta: " << stationarity << std::endl;
			std::cout << "constraints  eq: " << eq_con_values << std::endl;
			std::cout << "constraints ieq: " << ieq_con_values << std::endl;
		
			std::cout << "jac eq: " << nlp_->constraints_equality_jacobian(iterate) << std::endl;
		
			std::cout << "       dual  eq: " << iterate.dual_eq_values_ << std::endl;
			std::cout << "       dual ieq: " << iterate.dual_ieq_values_ << std::endl;
		}
		return resid_helper(iterate, stationarity, eq_con_values, ieq_con_values, iterate.dual_eq_values_, iterate.dual_ieq_values_);
	}
	double operator()(const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
		// first entry of rho(x,y,bary,mu):
		bool PRINT=true;
		std::vector<double> stationarity = nlp_->objective_gradient(iterate);
		
		matrix hessian = nlp_->lagrangian_hessian(iterate);
		
		std::vector<double> hessian_times_step = hessian.multiply(step.primal_values_);
		
		for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
			// mu*grad f + H*d:
			stationarity[stationarity_index] = iterate.penalty_parameter_*stationarity[stationarity_index] + hessian_times_step[stationarity_index];
			
		}
		
		matrix jacobian_eq = nlp_->constraints_equality_jacobian(iterate);
		std::vector<double> con_values_eq = nlp_->constraints_equality(iterate);
		// std::cout << "about to do j eq" << std::endl;
		std::vector<double> linear_step_eq = jacobian_eq.multiply(step.primal_values_);
		
		// std::cout << "about to do j ieq" << std::endl;
		matrix jacobian_ieq = nlp_->constraints_inequality_jacobian(iterate);
		std::vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
		std::vector<double> linear_step_ieq = jacobian_ieq.multiply(step.primal_values_);

		for (size_t eq_con_index=0; eq_con_index< nlp_->num_dual_eq(); ++eq_con_index) {
			con_values_eq[eq_con_index] += linear_step_eq[eq_con_index];
		}
		for (size_t ieq_con_index=0; ieq_con_index< nlp_->num_dual_ieq(); ++ieq_con_index) {
			con_values_ieq[ieq_con_index] += linear_step_ieq[ieq_con_index];
		}		
		return resid_helper(iterate, stationarity, con_values_eq, con_values_ieq, step.dual_eq_values_, step.dual_ieq_values_);
	}
protected:
};

class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(Nlp &nlp) : FunctionWithNLPState(nlp), example_(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual()), first_(true)  {}
	
	iSQOStep operator()(const iSQOQuadraticSubproblem &subproblem) {
		// qpOASES::SQProblem example_(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual());
		/* Setting up QProblem object. */
		int nWSR = 2000;
		qpOASES::returnValue ret;
		// std::cout << std::endl;
		// std::cout << "hessian: " << subproblem.hessian_ <<std::endl;
		// std::cout << "jacobian: " << subproblem.jacobian_ <<std::endl;
		// std::cout << "gradient: " << subproblem.gradient_ << std::endl;
		// std::cout << "lb: : " << subproblem.lower_bound_ << std::endl;
		// std::cout << "ub: : " << subproblem.upper_bound_ << std::endl;
		// std::cout << "lbA: " << subproblem.jacobian_lower_bound_ << std::endl;
		// std::cout << "ubA: " << subproblem.jacobian_upper_bound_ << std::endl;
		
		qpOASES::Options opt;
		// opt.setToReliable();
		// opt.enableRegularisation = qpOASES::BooleanType(true);

		opt.terminationTolerance = 1e-6;
		example_.setOptions(opt);
		example_.setPrintLevel(qpOASES::PL_NONE);
		// std::cout.flush();
		// std::cerr.flush();
		if (first_) {
			// std::cout << "subproblem.hessian_.data_[0]: " << subproblem.hessian_ << std::endl;
			// std::cout << "subproblem: " << subproblem.gradient_ << std::endl;
			// std::cout << "subproblem: " << subproblem.jacobian_ << std::endl;
			// std::cout << "subproblem: " << subproblem.lower_bound_ << std::endl;
			// std::cout << "subproblem: " << subproblem.upper_bound_ << std::endl;
			// std::cout << "subproblem: " << subproblem.jacobian_lower_bound_ << std::endl;
			// std::cout << "subproblem: " << subproblem.jacobian_upper_bound_ << std::endl;
			// subproblem.print();
			ret = example_.init( &subproblem.hessian_.data_[0],
								 &subproblem.gradient_[0],
								 &subproblem.jacobian_.data_[0],
								 &subproblem.lower_bound_[0],
								 &subproblem.upper_bound_[0],
								 &subproblem.jacobian_lower_bound_[0],
								 &subproblem.jacobian_upper_bound_[0],
								 nWSR,
								 0);
		} else{
			
			// subproblem.print();
			ret = example_.hotstart(&subproblem.hessian_.data_[0],
								 &subproblem.gradient_[0],
								 &subproblem.jacobian_.data_[0],
								 &subproblem.lower_bound_[0],
								 &subproblem.upper_bound_[0],
								 &subproblem.jacobian_lower_bound_[0],
								 &subproblem.jacobian_upper_bound_[0],
								 nWSR,
								 0);
		}
		std::cerr.flush();
		std::cout.flush();

		size_t return_status = example_.getStatus();
		// getGlobalMessageHandler()->listAllMessages();
		if( ret != qpOASES::SUCCESSFUL_RETURN ){
			// subproblem.print();
	        // printf( "%s\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );
			first_ = true;
			example_.reset();
		} else {
			if (first_) first_ = false;
		}
		iSQOStep step(subproblem.num_nlp_variables_,subproblem.num_nlp_constraints_eq_,subproblem.num_nlp_constraints_ieq_);
		// std::cout << "address of step: " << &step << std::endl;
		
		// full_primal needs to hold primal variables and all slack variables
		std::vector<double> primal_and_slack_values(nlp_->num_primal() + 2*nlp_->num_dual());
		example_.getPrimalSolution( &primal_and_slack_values[0] );
		for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index)
			step.primal_values_[primal_index] = primal_and_slack_values[primal_index];
		
		// std::cout << "primal: " << step.primal_values_ << std::endl;
		// 	- every QP variable has a multiplier:
		//		-- iterate.num_primal_: primal variables in NLP
		// 		-- 2*iterate.num_dual_eq_: positive & negative slacks for equalities
		//		-- 2*iterate.num_dual_ieq_: positive & negative slacks for inequalities
		//	- every constraint has a multiplier:
		//		-- iterate.num_dual_eq_
		//		-- iterate.num_dual_ieq_
		// so in total: 
		std::vector<double> yOpt(subproblem.num_qp_variables_ + subproblem.num_qp_constraints_);
		example_.getDualSolution( &yOpt[0] );				
		
		step.status_ = ret;
		// to get eq con multipliers, read past NLP & slack variables:
		for (size_t dual_eq_index=0; dual_eq_index<subproblem.num_nlp_constraints_eq_; ++dual_eq_index) {
			step.dual_eq_values_[dual_eq_index] = -yOpt[subproblem.num_qp_variables_+dual_eq_index];
		}
		// to get inequality constraint multipliers, read past NLP & slack variables and eq con multipliers:
		for (size_t dual_ieq_index=0; dual_ieq_index<subproblem.num_nlp_constraints_ieq_; ++dual_ieq_index) {
			step.dual_ieq_values_[dual_ieq_index] = -yOpt[subproblem.num_qp_variables_+subproblem.num_nlp_constraints_eq_+dual_ieq_index];
		}
		// std::cout << "dual eq: " << step.dual_eq_values_ << std::endl;
		// std::cout << "dual ieq: " << step.dual_ieq_values_ << std::endl;
		
		bool PRINT=false;
		if (PRINT) {
			std::cout << step;
		}
		return step;
	}
private:
protected:
	qpOASES::SQProblem example_;
	bool first_;
};

class LinearModelFunction : public FunctionWithNLPState {
public:
	LinearModelFunction(Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp) {
		
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) {
		double f = nlp_->objective(iterate);
		std::vector<double> gradient = nlp_->objective_gradient(iterate);
		// dotproduct(gradient, )
		double dot_product=0.0;
		for (size_t primal_index=0; primal_index<iterate.num_primal_; ++primal_index) {
			dot_product += gradient[primal_index]*step.primal_values_[primal_index];
		}
		return iterate.penalty_parameter_*(f + dot_product) + constraint_violation_func_(iterate, step);
	}
private:
protected:
	ConstraintViolationFunction constraint_violation_func_;
};

class LinearReductionFunction : public FunctionWithNLPState {
public:
	LinearReductionFunction(Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp) {
		
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) const {
		bool PRINT=false;
		// double f = nlp_->objective(iterate);
		std::vector<double> gradient = nlp_->objective_gradient(iterate);
		// dotproduct(gradient, )
		double dot_product=0.0;
		for (size_t primal_index=0; primal_index<iterate.num_primal_; ++primal_index) {
			dot_product += gradient[primal_index]*step.primal_values_[primal_index];
		}
		
		if (PRINT) std::cout << "linear decrease: " << std::endl;
		if (PRINT) std::cout << " - " << -iterate.penalty_parameter_*dot_product << std::endl;
		if (PRINT) std::cout << " - " << constraint_violation_func_(iterate) << std::endl;
		if (PRINT) std::cout << " - " << constraint_violation_func_(iterate,step) << std::endl;
		if (PRINT) std::cout << " - " << -iterate.penalty_parameter_*dot_product + constraint_violation_func_(iterate) - constraint_violation_func_(iterate,step) << std::endl;
		return -iterate.penalty_parameter_*dot_product + constraint_violation_func_(iterate) - constraint_violation_func_(iterate,step);
	}
private:
protected:
	ConstraintViolationFunction constraint_violation_func_;
};

class LineSearchFunction : public FunctionWithNLPState {
public:
	LineSearchFunction(Nlp &nlp) : FunctionWithNLPState(nlp), penalty_func_(nlp), linear_decrease_func_(nlp) {
		// std::cout << "initializing a linesearcher!" << std::endl;
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) {
		// bugs found here: it_and_step pulled from nlp->initial, which meant it didn't have the right penalty parameter, which caused the linesearch to fail.
		iSQOIterate it_and_step(iterate);
		bool PRINT = false;
		double step_size = 1.0;
		double penfunc_start = penalty_func_(iterate);
		double linear_reduction_start = linear_decrease_func_(iterate, step);
		double machine_precision = std::numeric_limits<double>::epsilon();
		
		double eta = 1e-8;
		double linesearch_step_size_reduction = 0.5;
		if (PRINT) std::cout << std::endl;
		for (size_t alpha_cuts = 0; alpha_cuts < 20; alpha_cuts++) {
			it_and_step.update(iterate, step_size, step);
			double penfunc_step = penalty_func_(it_and_step);
			
			if (PRINT) {
				std::cout << "step_size: " << step_size << std::endl;
				std::cout << "f: " << nlp_->objective(it_and_step) << std::endl;
				std::cout << "c: " << nlp_->constraints_equality(it_and_step) << std::endl;
				std::cout << "alpha_cut: " << alpha_cuts << ": original: " << penfunc_start << "; new: " << penfunc_step << "; reduction: " << (penfunc_start-penfunc_step) << std::endl;
				std::cout << "  rhs: " << - eta*step_size*linear_reduction_start  << " + " << 10*machine_precision*std::max(std::abs(penfunc_step), 1.0) << " = " << (-eta*step_size*linear_reduction_start+10*machine_precision*std::max(std::abs(penfunc_step), 1.0)) << std::endl;
				std::cout << "mu = " << it_and_step.penalty_parameter_ << std::endl;
			}
			if (penfunc_step - penfunc_start <= - eta*step_size*linear_reduction_start + 10*machine_precision*std::max(std::abs(penfunc_step), 1.0)) {
				break;
			}
			step_size = linesearch_step_size_reduction*step_size;
		}
		return step_size;
	}
private:
protected:
	PenaltyFunction penalty_func_;
	LinearReductionFunction linear_decrease_func_;
};
class HessianShifter : public FunctionWithNLPState {
public:
	HessianShifter(Nlp &nlp) : FunctionWithNLPState(nlp), solve_qp_(nlp), last_shift_(0.0) {
		
	}
	iSQOStep operator()(const iSQOIterate &iterate, iSQOQuadraticSubproblem &subproblem) {
		bool PRINT=false;
		double shift_w0 = 1e-4;
		double shift_min = 1e-20;
		double shiftkwbarp = 100;
		double shiftkwp = 8;
		double shift_kwm = 1.0/3.0;
		double shiftmax = 1e40;
		
		double current_shift = 0.0;
		
		// try to solve without regularization:
		subproblem.inc_regularization(current_shift);
		iSQOStep return_step = solve_qp_(subproblem);
		
		
		int current_regularization_steps = 0;
		std::vector<double> current_step_values(nlp_->num_primal() + 2*nlp_->num_dual());
		for (size_t primal_index = 0; primal_index < nlp_->num_primal(); ++primal_index) {
			current_step_values[primal_index] = return_step.primal_values_[primal_index];
		}
		
		// matrix hessian = nlp_->lagrangian_hessian(iterate);
		std::vector<double> hessian_step = subproblem.hessian_.multiply(current_step_values);
		if (PRINT) std::cout << "hessian step: " << subproblem.hessian_ << std::endl;
		double total = 0.0;
		for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index) {
			total += hessian_step[primal_index] * return_step.primal_values_[primal_index];
		}
		if (PRINT) std::cout << std::endl << "total: " << total << "; norm: " << 1e-8*return_step.x_norm()*return_step.x_norm() << std::endl;
		if (PRINT) std::cout << "required: " << (total < (1e-8*return_step.x_norm()*return_step.x_norm())) << std::endl;
		// std::vector<double> hessian_step = subproblem.hessian_.multiply()
		while ((return_step.status_ != 0) || (return_step.x_norm() > 1e9) || (total < (1e-8*return_step.x_norm()*return_step.x_norm()))) {
			if (current_regularization_steps == 0) {
				if (last_shift_ == 0.0) {
					current_shift = shift_w0;
				} else {
					current_shift = std::max(shift_min, shift_kwm * last_shift_);
				}
			} else {
				if (last_shift_ == 0.0) {
					current_shift = shiftkwbarp * current_shift;
				} else {
					current_shift = shiftkwp*current_shift;
				}
			}
			subproblem.inc_regularization(current_shift);
			return_step = solve_qp_(subproblem);
			
			for (size_t primal_index = 0; primal_index < nlp_->num_primal(); ++primal_index) {
				current_step_values[primal_index] = return_step.primal_values_[primal_index];
			}
			if (PRINT) std::cout << "hessian: " << subproblem.hessian_ << std::endl;
			if (PRINT) std::cout << "return step: " << return_step.primal_values_ << std::endl;
			hessian_step = subproblem.hessian_.multiply(current_step_values);
			if (PRINT) std::cout << "hessian step: " << hessian_step << std::endl;
			total = 0.0;
			for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index) {
				total += hessian_step[primal_index] * return_step.primal_values_[primal_index];
			}
			if (PRINT) std::cout << "current_shift: " << current_shift << std::endl;
			if (PRINT) std::cout << std::endl << "total: " << total << "; norm: " << 1e-8*return_step.x_norm()*return_step.x_norm() << std::endl;
			if (PRINT) std::cout << "required: " << (total < (1e-8*return_step.x_norm()*return_step.x_norm())) << std::endl;
			
			++current_regularization_steps;
			
			if (current_regularization_steps == 20){
				std::cout << "FAILURE!" << std::endl;
				assert(false);
				break;
			}
		}
		last_shift_ = current_shift;
		return return_step;
	}
	
	double get_last_shift() const {
		return last_shift_;
	}
private:
protected:
	SolveQuadraticProgram solve_qp_;
	double last_shift_;
};

class TextOutput : public FunctionWithNLPState {
public:
	TextOutput (Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp),linear_decrease_func_(nlp), pen_func_(nlp), residual_func_(nlp) {
		
	}
	
	void start() {
		printf("-----|");
		printf("----------------------|");
		printf("----------------------|");
		printf("----------------------&");
		printf("-----------------[ FEASIBILITY ]-----------------|");
		printf("-------------------[ PENALTY ]-------------------|");
		printf("--------------------------------------------|");
		printf("----------\n");
		printf("%s",output_desc_pre_);
		printf("%s",output_desc_subprob_);
		printf("%s",output_desc_subprob_);
		printf("%s",output_desc_post_);
		printf("-----|");
		printf("----------------------|");
		printf("----------------------|");
		printf("----------------------&");
		printf("-------------------------------------------------|");
		printf("-------------------------------------------------|");
		printf("--------------------------------------------|");
		printf("----------\n");
	}
	void pre(size_t iter, const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate) const {
		// iter, problem, feas iterate, pen iterate, constraintviolation, iterate, penfunc, residual_func
		//  -> params: iter, feas iter, pen iter
		//	-> state: nlp_, constraint_violation_func_, pen_func_, residual_func_
		printf(output_format_pre_, 
				iter, 
				nlp_->objective(penalty_iterate), constraint_violation_func_(penalty_iterate),
				penalty_iterate.penalty_parameter_, pen_func_(penalty_iterate),
				residual_func_(feasibility_iterate), residual_func_(penalty_iterate)
				);
	}
	void subproblem(double shift, const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
		// shift, step, linear_decrease_func, residual_func
		//  -> params: shift, step
		//	-> state: linear_decrease_func, residual_func_
		printf(output_format_subprob_,
				shift, step.status_, step.x_norm(),
				linear_decrease_func_(iterate, step), residual_func_(iterate, subproblem, step)
					);
	}
	void subproblem_skip() {
		printf("         -   -           -          -          - |");
	}
	void post(const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate, const iSQOStep &combination_step, std::string step_type, double step_mix, double alpha) {
		// step, feas_iter, penalty_iter, alpha
		//  -> params: step, feas_iter, penalty_iter, alpha
		printf(output_format_post_,
				step_type.c_str(),
				step_mix, combination_step.x_norm(), linear_decrease_func_(feasibility_iterate,combination_step), linear_decrease_func_(penalty_iterate, combination_step),
				alpha
					);
	}
protected:
	// Nlp *nlp_;
	ConstraintViolationFunction constraint_violation_func_;
	LinearReductionFunction linear_decrease_func_;
	PenaltyFunction pen_func_;
	ResidualFunction residual_func_;
	
	static const char output_desc_pre_[];
	static const char output_desc_subprob_[];
	static const char output_desc_post_[];
	static const char output_format_pre_[];
	static const char output_format_subprob_[];
	static const char output_format_post_[];
};
const char TextOutput::output_desc_pre_[] = " it  |       obj     infeas |       pen      merit |   feaskkt     penkkt &";
const char TextOutput::output_desc_subprob_[] = "     shift  msg     ||d||     penred        res  |";
const char TextOutput::output_desc_post_[] = " TT   CvxComb    ||d||   FeasRed     PenRed |    alpha\n";
const char TextOutput::output_format_pre_[] = " %3d | %9.2e  %9.2e | %9.2e  %+9.2e | %9.2e  %9.2e &";
const char TextOutput::output_format_subprob_[] = " %9.2e  %3d  %9.2e  %9.2e  %9.2e |";
const char TextOutput::output_format_post_[] = " %s %9.2e %9.2e %9.2e %9.2e |%9.2e\n";


bool assert_close(double val1, double val2, double tol) {
	assert(abs(val1) <= abs(val2) + tol*abs(val1));
	assert(abs(val2) <= abs(val1) + tol*abs(val2));
	assert(((val1 > tol) && (val2 > tol)) || ((val1 < -tol) && (val2 < -tol)) || (val1 == 0 && val2 == 0));
	return (abs(val1) <= abs(val2) + tol*abs(val1)) && (abs(val2) <= abs(val1) + tol*abs(val2));
}
int main(int argc, char **argv) {
	std::string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/hs067.nl");
	if (argc>1) {
		problem_file = std::string(argv[1]);
	}
	AmplNlp problem(problem_file);
	iSQOIterate penalty_iterate = problem.initial();
    // penalty_iterate.primal_values_[0] = 1.25*3;
    // penalty_iterate.primal_values_[1] = 1.25*3;
    // penalty_iterate.primal_values_[2] = 1.25*3;
    // penalty_iterate.primal_values_[3] = 1.25*3;
    // penalty_iterate.primal_values_[4] = 1.25*3;
    // penalty_iterate.primal_values_[5] = 1.25*3;
    // penalty_iterate.primal_values_[7] = 1.25*3;
    // penalty_iterate.primal_values_[8] = 1.25*3;
    // penalty_iterate.primal_values_[9] = 1.25*3;
    
	
    std::cout << "objective: " << problem.objective(penalty_iterate) << std::endl;
	std::cout << "constraints: "<< problem.constraints_equality(penalty_iterate) << std::endl;
	
    qpOASES::SparseMatrix sparse_jacobian = problem.constraints_equality_jacobian_sparse(penalty_iterate);
    // std::cout << "eq jacobian: " << ) << std::endl;
    
}
int main3() {
	std::string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/hs009.nl");
	AmplNlp problem(problem_file);
	iSQOIterate penalty_iterate = problem.initial();
	matrix je = problem.constraints_equality_jacobian(penalty_iterate);
	std::vector<double> x(penalty_iterate.primal_values_);
	std::vector<double> y(penalty_iterate.dual_eq_values_);
	std::cout << "jac: "<< je << std::endl;
	std::cout << "x: " << x << std::endl;
	std::cout << je.multiply(x) << std::endl;
	
	std::cout << "y: " << y << std::endl;
	std::cout << je.multiply_transpose(y) << std::endl;
	
}
int main2() {
	std::cout << "okay, so..." << std::endl;
	
	// tests for constraint conversation:
	// equality constrained problems: (to c(x) = 0)
	//	* 0 <= c_1(x) <= 0  ---> A_1d = 0 - c_1(x_k)
	//	* 1 <= c_2(x) <= 1  ---> A_2d = 1 - c_2(x_k)
	// inequality constrained problems: (to barc(x) <= 0)
	//  * double constrained: 	-1 <= c_3(x) <= 5  --> A_3d <= 5 - c_3(x_k) AND -A_3d <= -(-1) - (-c_3(x_k))
	//			c_3(x) <= 5  ==> 
	//			-1 <= c_3(x) ==> -c_3(x) <= -(-1) ==> -A_3d - c_3(x_k) <= -(-1) ==> -A_3d  <= -(-1) - (-c_3(x_k))
	//	* double constrained: 	0 <= c_4(x) <= 5 --> A_4d <= 5 - c_4(x_k) AND -A_4d <= -(0) - (-c_4(x_k))
	//	* double constrained: 	-1 <= c_5(x) <= 0 --> A_5d <= 0 - c_5(x_k) AND -A_5d <= -(-1) - (-c_5(x_k))
	//	* only upper: -INF <= c_6(x) <= 0 --> A_6d <= 0 - c_6(x_k)
	//	* only upper: -INF <= c_7(x) <= 5 --> A_7d <= 5 - c_7(x_k)
	//	* only lower: 0 <= c_8(x) <= INF --> -A_8d <= -(0) - c_8(x_k)
	//	* only lower: 5 <= c_9(x) <= INF --> -A_9d <= -(5) - c_9(x_k)
	
	// 
	// tests for variable conversion:
	// 
	// AmplNlp ampl_nlp_object("/Users/traviscj/optimization/cute_nl_nopresolve/hs016.nl");
	AmplNlp ampl_nlp_object("/Users/traviscj/ampl_bound_reformulation.nl");
	iSQOIterate iterate = ampl_nlp_object.initial();
	std::cout << "initial iterate: " << iterate.primal_values_ << std::endl;
	// std::cout << "Zeroth Order:" << std::endl;
	std::cout << "===[Objective]==============" << std::endl;
	double ampl_obj = ampl_nlp_object.objective(iterate);
	assert(ampl_obj == 5.921590000000000e+00);//, "ampl_obj wrong");
	std::cout << "  AMPL : " << ampl_nlp_object.objective(iterate) << std::endl;
	
	std::cout << "===[Objective Gradient]==============" << std::endl;
	std::vector<double> ampl_obj_gradient = ampl_nlp_object.objective_gradient(iterate);
	assert(ampl_obj_gradient[0] == 1.0e0);
	assert(ampl_obj_gradient[1] == 1.0e0);
	std::cout << "  AMPL : " << ampl_obj_gradient << std::endl;
	
	
	std::cout << "===[Equality Constraints]==============" << std::endl;
	std::vector<double> con_eq_values = ampl_nlp_object.constraints_equality(iterate);
	std::cout << "  AMPL : " << con_eq_values << std::endl;
	assert(con_eq_values[0]==-41.924375456199996); // negative because switches var to left-hand-side.
	assert(con_eq_values[1]==141.09951409669998);
	
	std::cout << "===[Equality Jacobian]==============" << std::endl;
	matrix eq_jacobian = ampl_nlp_object.constraints_equality_jacobian(iterate);
	std::cout << "  AMPL : " << eq_jacobian << std::endl;
	assert(eq_jacobian.get(0,0) == -1.256636000000000e+01);
	assert(eq_jacobian.get(0,1) == -1.668000000000000e+01);
	assert(eq_jacobian.get(1,0) == 4.398226000000000e+01);
	assert(eq_jacobian.get(1,1) == 6.116000000000000e+01);
	
	
	
	
	std::cout << "===[Inequality Constraints]==============" << std::endl;
	std::vector<double> ieq_values = ampl_nlp_object.constraints_inequality(iterate);
	std::cout << "  AMPL : " << ieq_values << std::endl;
	double ieq_tol = 1e-12;
	assert(abs(ieq_values[0]) <= (abs(-3.482753668339000e+02) + ieq_tol*abs(ieq_values[0])));
	assert(abs(ieq_values[1]) <= (abs(-6.820391459396999e+02) + ieq_tol*abs(ieq_values[1])));
	assert(abs(ieq_values[2]) <= (abs(3.362753668339000e+02) + ieq_tol*abs(ieq_values[2])));
	assert(abs(ieq_values[3]) <= (abs(7.346270723082998e+02) + ieq_tol*abs(ieq_values[3])));
	
	// -6.820391459396999e+02
	
	std::cout << "===[Inequality Jacobian]==============" << std::endl;
	matrix ieq_jacobian = ampl_nlp_object.constraints_inequality_jacobian(iterate);
	std::cout << "  AMPL : " << ieq_jacobian << std::endl;
	assert_close(ieq_jacobian.get(0,0), -1.193804200000000e+02, 1e-10); // c2 lower
	assert_close(ieq_jacobian.get(0,1), -1.278800000000000e+02, 1e-10);
	assert_close(ieq_jacobian.get(1,0), -2.324776600000000e+02, 1e-10); // c3 lower
	assert_close(ieq_jacobian.get(1,1), -2.279600000000000e+02, 1e-10);
	assert_close(ieq_jacobian.get(2,0), 1.193804200000000e+02, 1e-10);  // c2 upper
	assert_close(ieq_jacobian.get(2,1), 1.278800000000000e+02, 1e-10);  
	assert_close(ieq_jacobian.get(3,0), 2.701767400000000e+02, 1e-10);  // c4 upper
	assert_close(ieq_jacobian.get(3,1), 2.613200000000000e+02, 1e-10);
	
	
	
	
	
	
	// assert(ieq_jacobian.get(3,0) == 1.193804200000000e+02); assert(ieq_jacobian.get(3,1) == 1.278800000000000e+02);  // c2 upper
	
	std::cout << std::endl;
	std::cout << "TESTED FRONTIER!" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	
	std::cout << "===[Lagrangian Hessian]==============" << std::endl;
	matrix lh_f_only = ampl_nlp_object.lagrangian_hessian(iterate);
	std::cout << "  AMPL : " << lh_f_only << std::endl;
	assert_close(lh_f_only.get(0,0), 0.00e0, 1e-10);  // c2 upper
	assert_close(lh_f_only.get(0,1), 0.00e0, 1e-10);  
	assert_close(lh_f_only.get(1,0), 0.00e0, 1e-10);  // c4 upper
	assert_close(lh_f_only.get(1,1), 0.00e0, 1e-10);
	
	iterate.dual_eq_values_[0] = 1.0;
	matrix lh_c0_only = ampl_nlp_object.lagrangian_hessian(iterate);
	std::cout << "  AMPL : " << lh_c0_only << std::endl;
	iterate.dual_eq_values_[0] = -1.0;
	matrix lh_c0_only_flip = ampl_nlp_object.lagrangian_hessian(iterate);
	std::cout << "  AMPL : " << lh_c0_only_flip << std::endl;
	
	
	return 0;
}
int main1() {
	AmplNlp ampl_nlp_object("/Users/traviscj/optimization/cute_nl_nopresolve/hs014.nl");
	
	
	Hs014 problem;
	iSQOIterate iterate(2,1,1);
	iterate.penalty_parameter_ = 1.0e-1;
	iterate.primal_values_[0] = 2.0;
	iterate.primal_values_[1] = 2.0;
	iterate.dual_eq_values_[0] = -1.0;
	iterate.dual_ieq_values_[0] = 1000;
	// iterate.dual_ieq_values_[0] = ;
	// iterate.dual_ieq_values_[0] = 0;
	// iterate.penalty_parameter_ = 1.0e-1;
	// iterate.primal_values_[0] = .822875656;
	// iterate.primal_values_[1] = .911437828;
	// iterate.dual_eq_values_[0] = 1.846589027861980e+00;
	// iterate.dual_ieq_values_[0] = 1.594493103554523e+00;

	std::cout << "Zeroth Order:" << std::endl;
	std::cout << "===[Objective]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.objective(iterate) << std::endl;
	std::cout << "Custom : " << problem.objective(iterate) << std::endl;
	
	std::cout << "===[Equality Constraints]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.constraints_equality(iterate) << std::endl;
	std::cout << "Custom : " << problem.constraints_equality(iterate) << std::endl;
	
	std::cout << "===[Inequality Constraints]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.constraints_inequality(iterate) << std::endl;
	std::cout << "Custom : " << problem.constraints_inequality(iterate) << std::endl;
	
	std::cout << "First Order:" << std::endl;
	std::cout << "===[Objective Gradient]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.objective_gradient(iterate) << std::endl;
	std::cout << "Custom : " << problem.objective_gradient(iterate) << std::endl;
	
	std::cout << "===[Equality Jacobian]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.constraints_equality_jacobian(iterate) << std::endl;
	std::cout << "Custom : " << problem.constraints_equality_jacobian(iterate) << std::endl;
	
	
	std::cout << "===[Inequality Jacobian]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.constraints_inequality_jacobian(iterate) << std::endl;
	std::cout << "Custom : " << problem.constraints_inequality_jacobian(iterate) << std::endl;
	
	std::cout << "Second Order:" << std::endl;
	std::cout << "===[Lagrangian Hessian]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.lagrangian_hessian(iterate) << std::endl;
	std::cout << "Custom : " << problem.lagrangian_hessian(iterate) << std::endl;
		
	return 0;
}
int main0(int argc, char **argv) {
	
	std::string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/hs067.nl");
	if (argc>1) {
		problem_file = std::string(argv[1]);
	}
	AmplNlp problem(problem_file);
	
	// Utilities for NLP:
	PenaltyFunction penalty_function(problem);
	ResidualFunction kkt_residual(problem);
	SolveQuadraticProgram solve_qp(problem);
	ConstraintViolationFunction constraint_violation(problem);
	LineSearchFunction line_search(problem);
	LinearReductionFunction linear_reduction(problem);
	HessianShifter hessian_shifting_penalty_qp_solve(problem);
	HessianShifter hessian_shifting_feasibility_qp_solve(problem);
	// Other
	TextOutput text_output(problem);
	
	// Paper-defined Parameter Values
	double linear_decrease_threshold = 1e-1; // epsilon in the paper.
	double linear_reduction_threshold_for_penalty_reduction = 1e-2; // beta in the paper.
	double mostly_feasibility_step_threshold = 1e-3; // tau in the paper.
	double penalty_parameter_reduction_factor = 0.2; // delta in paper.
	double termination_threshold = 1e-6;
	double termination_penalty_threshold = 1e-8;
	double hessian_min_convexity = 1e-8;	// theta in the paper.
	// Code-defined parameter values:
	double convex_combination_search_decrease = .9;
	bool PRINT=false;
	size_t maximum_iterations = 1000;
	size_t max_num_comb_reductions = 100;
	double machine_precision = std::numeric_limits<double>::epsilon();
	
	text_output.start();
	
	//////////////////////////
	// ALGORITHM A // STEP 1
	//////////////////////////
	iSQOIterate penalty_iterate = problem.initial();
	penalty_iterate.penalty_parameter_ = 1.0e-1;

	iSQOIterate feasibility_iterate = problem.initial();
	feasibility_iterate.penalty_parameter_ = 0.0;
	size_t iter=-1;
	for (iter = 0; iter < maximum_iterations; iter ++ ) {
		iSQOStep combination_step(problem.num_primal(), problem.num_dual_eq(), problem.num_dual_ieq());
		
		text_output.pre(iter, feasibility_iterate, penalty_iterate);
		
		//////////////////////////
		// ALGORITHM A // STEP 2
		//////////////////////////
		if ((kkt_residual(penalty_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) < termination_threshold)) {
			std::cout << std::endl << "Termination 2a - optimality" << std::endl;
			break;
		}
		if ((kkt_residual(feasibility_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) > termination_threshold) && (penalty_iterate.penalty_parameter_ < termination_penalty_threshold)) {
			std::cout << std::endl << "Termination 2b - infeasible problem!" << std::endl;
			break;
		}
		
		//////////////////////////
		// ALGORITHM A // STEP 3
		//////////////////////////
		// Penalty problem is set up AND SOLVED, 
		iSQOQuadraticSubproblem penalty_subproblem(problem, penalty_iterate);
		iSQOStep penalty_step = hessian_shifting_penalty_qp_solve(penalty_iterate, penalty_subproblem);
		iSQOStep feasibility_step(problem.num_primal(),problem.num_dual_eq(),problem.num_dual_ieq());
		
		// per-iteration variable setup.
		std::string steptype = "4a";
		double combination_step_contribution_from_penalty_step = -1.0;
		
		if (linear_reduction(penalty_iterate,penalty_step) >= linear_decrease_threshold*constraint_violation(penalty_iterate) + 10*machine_precision*linear_reduction(penalty_iterate,penalty_step)) {
			//////////////////////////
			// ALGORITHM A // STEP 3a
			combination_step_contribution_from_penalty_step = 1.0;
			text_output.subproblem_skip();

			combination_step.set_primal(penalty_step);
		} else {
			
			//////////////////////////
			// ALGORITHM A // STEP 4
			//////////////////////////
			// Feasibility problem is set up and solved.
			iSQOQuadraticSubproblem feasibility_subproblem(problem, feasibility_iterate);
			feasibility_step = hessian_shifting_feasibility_qp_solve(feasibility_iterate, feasibility_subproblem);
			
			if (linear_reduction(penalty_iterate,penalty_step) >= linear_decrease_threshold*linear_reduction(feasibility_iterate,feasibility_step)) {
				//////////////////////////
				// ALGORITHM A // STEP 4a
				combination_step_contribution_from_penalty_step = 1.0;
				steptype = "5a";
			} else {
				//////////////////////////
				// ALGORITHM A // STEP 4b
				combination_step_contribution_from_penalty_step = 1.0;
				
				int num_comb_reductions = 0;
				combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
				while (! (linear_reduction(feasibility_iterate,combination_step)  >=linear_decrease_threshold*linear_reduction(feasibility_iterate,feasibility_step) + 10*machine_precision*linear_reduction(feasibility_iterate,combination_step)) ) {
					combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
					
					combination_step_contribution_from_penalty_step = convex_combination_search_decrease*combination_step_contribution_from_penalty_step;
					++num_comb_reductions;
					assert(! (num_comb_reductions >= max_num_comb_reductions));
				}
				steptype = "5b";
				
				// beta: 
				double linear_decrease_in_penalty_combination = linear_reduction(penalty_iterate,combination_step);
				double linear_decrease_in_feasibility_combination = linear_reduction(feasibility_iterate,combination_step);
				if (linear_decrease_in_penalty_combination >= linear_reduction_threshold_for_penalty_reduction*linear_decrease_in_feasibility_combination + 10*machine_precision*linear_decrease_in_penalty_combination) {
					if (combination_step_contribution_from_penalty_step >= mostly_feasibility_step_threshold) {
						// no-op: keep the current penalty parameter.
					} else {
						penalty_iterate.penalty_parameter_ = penalty_parameter_reduction_factor*penalty_iterate.penalty_parameter_;
					}
				} else {
					std::vector<double> gradient = problem.objective_gradient(penalty_iterate);
					double potential_new_penalty_parameter_numerator = (1-linear_reduction_threshold_for_penalty_reduction)*linear_decrease_in_feasibility_combination;
					double potential_new_penalty_parameter_denominator = combination_step.x_dot_product(gradient) + hessian_min_convexity*pow(combination_step.x_norm(),2);
					penalty_iterate.penalty_parameter_ = std::min(penalty_parameter_reduction_factor*penalty_iterate.penalty_parameter_, potential_new_penalty_parameter_numerator / potential_new_penalty_parameter_denominator);
				}
			}

			// combination gets the convex combination of penalty and feasibility:
			combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
			
			text_output.subproblem(hessian_shifting_feasibility_qp_solve.get_last_shift(), feasibility_iterate, feasibility_subproblem, feasibility_step);
		}
		
		//////////////////////////
		// ALGORITHM A // STEP 5
		//////////////////////////
		double step_size = line_search(penalty_iterate, combination_step);
		
		text_output.subproblem(hessian_shifting_penalty_qp_solve.get_last_shift(), penalty_iterate, penalty_subproblem, penalty_step);
		text_output.post(feasibility_iterate, penalty_iterate, combination_step, steptype, combination_step_contribution_from_penalty_step, step_size);
		
		//////////////////////////
		// ALGORITHM A // STEP 6
		//////////////////////////
		// update penalty iterate & dual values:
		penalty_iterate.update(penalty_iterate, step_size, combination_step);
		penalty_iterate.update_dual(penalty_step);
		// update feasibility iterate & dual values:
		feasibility_iterate.update(penalty_iterate, step_size, combination_step);
		feasibility_iterate.update_dual(feasibility_step);
	}
	
	if (iter == maximum_iterations) {
		std::cout << "Failure - Did not converge" << std::endl;
	}
	std::cout << penalty_iterate.primal_values_ << std::endl;
	return 0;
}


