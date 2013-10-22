// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)

#include <cstdio>
#include "text_output.hh"

const char TextOutput::output_desc_pre_[] = " it  |       obj     infeas |       pen      merit |   feaskkt     penkkt &";
const char TextOutput::output_desc_subprob_[] = "     shift  msg     ||d||     penred        res  |";
const char TextOutput::output_desc_post_[] = " TT   CvxComb    ||d||   FeasRed     PenRed |    alpha\n";
const char TextOutput::output_format_pre_[] = " %3d | %9.2e  %9.2e | %9.2e  %+9.2e | %9.2e  %9.2e &";
const char TextOutput::output_format_subprob_[] = " %9.2e  %3d  %9.2e  %9.2e  %9.2e |";
const char TextOutput::output_format_post_[] = " %s %9.2e %9.2e %9.2e %9.2e |%9.2e\n";

TextOutput::TextOutput (Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp),linear_decrease_func_(nlp), pen_func_(nlp), residual_func_(nlp) {
	
}

void TextOutput::nlp() {
    printf("nvar: %5d\tncon: %5d\n", nlp_->num_primal(), nlp_->num_dual());
}
void TextOutput::start() {
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
void TextOutput::pre(size_t iter, const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate) const {
	// iter, problem, feas iterate, pen iterate, constraintviolation, iterate, penfunc, residual_func
	//  -> params: iter, feas iter, pen iter
	//	-> state: nlp_, constraint_violation_func_, pen_func_, residual_func_
	printf(output_format_pre_, 
			iter, 
			nlp_->objective(penalty_iterate), constraint_violation_func_(penalty_iterate),
			penalty_iterate.get_penalty_parameter(), pen_func_(penalty_iterate),
			residual_func_(feasibility_iterate), residual_func_(penalty_iterate)
			);
}
void TextOutput::subproblem(double shift, const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
	// shift, step, linear_decrease_func, residual_func
	//  -> params: shift, step
	//	-> state: linear_decrease_func, residual_func_
	printf(output_format_subprob_,
			shift, step.get_status(), step.x_norm(),
			linear_decrease_func_(iterate, step), residual_func_(iterate, subproblem, step)
				);
}
void TextOutput::subproblem_skip() {
	printf("         -   -           -          -          - |");
}
void TextOutput::post(const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate, const iSQOStep &combination_step, std::string step_type, double step_mix, double alpha) {
	// step, feas_iter, penalty_iter, alpha
	//  -> params: step, feas_iter, penalty_iter, alpha
	printf(output_format_post_,
			step_type.c_str(),
			step_mix, combination_step.x_norm(), linear_decrease_func_(feasibility_iterate,combination_step), linear_decrease_func_(penalty_iterate, combination_step),
			alpha
				);
}