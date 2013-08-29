#ifndef HESSIAN_SHIFTER_HH_Q33I3P8W
#define HESSIAN_SHIFTER_HH_Q33I3P8W

class HessianShifter : public FunctionWithNLPState {
public:
	HessianShifter(Nlp &nlp);
	iSQOStep operator()(iSQOQuadraticSubproblem &subproblem);
	iSQOStep operator()(iSQOSparseQuadraticSubproblem &subproblem);
	double get_last_shift() const ;
private:
protected:
	SolveQuadraticProgram solve_qp_;
	double last_shift_;
};


#endif /* end of include guard: HESSIAN_SHIFTER_HH_Q33I3P8W */
