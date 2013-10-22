// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)
#ifndef HESSIAN_SHIFTER_HH_Q33I3P8W
#define HESSIAN_SHIFTER_HH_Q33I3P8W

#include <string>
#include <sstream>

#include "linear_model_reduction.hh"

//! \brief a basic hessian-shifting routine
//!
//! If \f$ \frac{1}{2} d^T H d < \theta \| d \|_2^2\f$, 
//! then \f$ H \gets H + \xi I \f$ and try to re-solve the QP.
class HessianShifter : public FunctionWithNLPState {
public:
	HessianShifter(iSQOControlPanel &control, Nlp &nlp);
    //! \brief solve the QP with appropriate shifts until we get a \f$d\f$ matching the above condition.
	iSQOStep operator()(const iSQOIterate &iterate, iSQOQuadraticSubproblem &subproblem);
    
    //! \brief allow clients to find the last shift.
	double get_last_shift() const ;
    
    void save_qp_state() { solve_qp_.save_qp_state(); }
    void restore_qp_state() { solve_qp_.restore_qp_state(); }
    
    std::string get_info_str() { return shifter_info.str(); }
private:
protected:
	SolveQuadraticProgram solve_qp_;
    LinearReductionFunction linear_model_reduction_;
	double last_shift_;
    std::stringstream shifter_info;
};


#endif /* end of include guard: HESSIAN_SHIFTER_HH_Q33I3P8W */
