// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)
#ifndef HESSIAN_SHIFTER_HH_Q33I3P8W
#define HESSIAN_SHIFTER_HH_Q33I3P8W

//! \brief a basic hessian-shifting routine
//!
//! If \f$ \frac{1}{2} d^T H d < \theta \| d \|_2^2\f$, 
//! then \f$ H \gets H + \xi I \f$ and try to re-solve the QP.
class HessianShifter : public FunctionWithNLPState {
public:
	HessianShifter(Nlp &nlp);
    //! \brief solve the QP with appropriate shifts until we get a \f$d\f$ matching the above condition.
	iSQOStep operator()(iSQOQuadraticSubproblem &subproblem);
    
    //! \brief allow clients to find the last shift.
	double get_last_shift() const ;
private:
protected:
	SolveQuadraticProgram solve_qp_;
	double last_shift_;
};


#endif /* end of include guard: HESSIAN_SHIFTER_HH_Q33I3P8W */
