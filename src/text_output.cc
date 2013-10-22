// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)

#include <cstdio>
#include "text_output.hh"

const char TextOutput::output_desc_pre_[] = " it  |       obj     infeas |       pen      merit |   feaskkt     penkkt &";
const char TextOutput::output_desc_subprob_[] = "     shift  msg  pivots   ||d||     penred        res  |";
const char TextOutput::output_desc_post_[] = " TT   CvxComb    ||d||   FeasRed     PenRed |";
const char TextOutput::output_desc_line_search_[] = "    alpha\n";
const char TextOutput::output_format_pre_[] = " %3d | %9.2e  %9.2e | %9.2e  %+9.2e | %9.2e  %9.2e &";
const char TextOutput::output_format_subprob_[] = " %9.2e  %3d  %5d %9.2e  %9.2e  %9.2e |";
const char TextOutput::output_format_post_[] = " %s %9.2e %9.2e %9.2e %9.2e |";
const char TextOutput::output_format_line_search_[] = "%9.2e\n";

TextOutput::TextOutput (iSQOControlPanel &control, Nlp &nlp) : FunctionWithNLPState(control, nlp), constraint_violation_func_(control, nlp),linear_decrease_func_(control, nlp), pen_func_(control, nlp), residual_func_(control, nlp) {
	text_output_file_ = fopen("my_output_file", "w");
}
TextOutput::~TextOutput() {
    fclose(text_output_file_);
}

void TextOutput::nlp() {
    fprintf(text_output_file_, "nvar: %5d\tncon: %5d (eq: %5d, ieq: %5d)\n", nlp_->num_primal(), nlp_->num_dual(), nlp_->num_dual_eq(), nlp_->num_dual_ieq());
}

void TextOutput::start() {
	fprintf(text_output_file_, "-----|");
	fprintf(text_output_file_, "----------------------|");
	fprintf(text_output_file_, "----------------------|");
	fprintf(text_output_file_, "----------------------&");
	fprintf(text_output_file_, "--------------------[ FEASIBILITY ]--------------------|");
	fprintf(text_output_file_, "----------------------[ PENALTY ]----------------------|");
	fprintf(text_output_file_, "--------------------------------------------|");
	fprintf(text_output_file_, "----------\n");
	fprintf(text_output_file_, "%s",output_desc_pre_);
	fprintf(text_output_file_, "%s",output_desc_subprob_);
	fprintf(text_output_file_, "%s",output_desc_subprob_);
	fprintf(text_output_file_, "%s",output_desc_post_);
    fprintf(text_output_file_, "%s",output_desc_line_search_);
	fprintf(text_output_file_, "-----|");
	fprintf(text_output_file_, "----------------------|");
	fprintf(text_output_file_, "----------------------|");
	fprintf(text_output_file_, "----------------------&");
	fprintf(text_output_file_, "-------------------------------------------------------|");
	fprintf(text_output_file_, "-------------------------------------------------------|");
	fprintf(text_output_file_, "--------------------------------------------|");
	fprintf(text_output_file_, "----------\n");
}
void TextOutput::pre(size_t iter, const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate) const {
	// iter, problem, feas iterate, pen iterate, constraintviolation, iterate, penfunc, residual_func
	//  -> params: iter, feas iter, pen iter
	//	-> state: nlp_, constraint_violation_func_, pen_func_, residual_func_
	fprintf(text_output_file_, output_format_pre_, 
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
	fprintf(text_output_file_, output_format_subprob_,
			shift, step.get_status(), step.get_pivots(), step.x_norm(),
			linear_decrease_func_(iterate, step), residual_func_(iterate, subproblem, step)
				);
}
void TextOutput::subproblem_skip() {
	fprintf(text_output_file_, "         -    -      -         -          -          - |");
}
// text_output.post(feasibility_iterate, penalty_iterate, combination_step, steptype, combination_step_contribution_from_penalty_step);

void TextOutput::post(const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate, const iSQOStep &combination_step, std::string step_type, double step_mix) {
	// step, feas_iter, penalty_iter, alpha
	//  -> params: step, feas_iter, penalty_iter, alpha
    std::cout << "**** computing feas decrease: " << std::endl;
    double feas_dec = linear_decrease_func_(feasibility_iterate,combination_step);
    std::cout << "**** computing feas decrease: " << std::endl;
    double pen_dec = linear_decrease_func_(penalty_iterate, combination_step);
    std::cout << "**** computing feas decrease: " << std::endl;
	fprintf(text_output_file_, output_format_post_,
			step_type.c_str(),
			step_mix, combination_step.x_norm(), feas_dec, pen_dec
				);
}

void TextOutput::line_search(double alpha) {
	// step, feas_iter, penalty_iter, alpha
	//  -> params: step, feas_iter, penalty_iter, alpha
	fprintf(text_output_file_, output_format_line_search_,
			alpha
				);
}