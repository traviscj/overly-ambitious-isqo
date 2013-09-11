
#include <iostream>
#include <memory>
#include <cassert>

#include "nlp_ampl.hh"
#include "matrix.hh"

AmplNlp::AmplNlp(std::string stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
    // , eq_jacobian_(0), ieq_jacobian_(0), hessian_(0)
    ConstructHelper(stub_str);
}

DenseAmplNlp::DenseAmplNlp(std::string stub_str) : AmplNlp(stub_str) {
    // ConstructHelper(stub_str);
}
DenseAmplNlp::~DenseAmplNlp() {}
SparseAmplNlp::SparseAmplNlp(std::string stub_str) : AmplNlp(stub_str) {
    // ConstructHelper(stub_str);
}

SparseAmplNlp::~SparseAmplNlp() {}

// AmplNlp::AmplNlp(std::string stub_str) : Nlp(-1,-1,-1), PRINT_(false), eq_jacobian_(0,0,0), ieq_jacobian_(0,0,0), hessian_(0,0,0) {
//     ConstructHelper(stub_str);
// }

// AmplNlp(char *stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
void AmplNlp::ConstructHelper(std::string stub_str) {
	if (PRINT_) std::cout << "Constructing an AmplNlp" << std::endl;
	ASL *asl;
	asl = ASL_alloc(ASL_read_pfgh);
	nerror_ = new int;
	*nerror_ = 0;
	// char * punt = stub_str.data();
	// if (PRINT_) 
		std::cout << "Filename: " << stub_str << std::endl;
	nl_ = jac0dim(&stub_str[0], (stub_str).length());
	if (PRINT_) 
		std::cout << "ran jac0dim: " << nl_ << std::endl;
		
	this->num_primal_ = n_var;
	X0_.resize(this->num_primal());
	X0 = &X0_[0];
	// A_vals = new double[nzc];
	// A_rownos = new int[nzc];
	// A_colstarts = new int[nzc];
	// ideally, we would use these....

	int status = pfgh_read(nl_, ASL_return_read_err | ASL_findgroups);

    assert(n_obj == 1);

	if (PRINT_) std::cout << "X0: " << X0[0] << ", " << X0[1] << std::endl;
	
	for (size_t ampl_constraint_index=0; ampl_constraint_index<n_con; ++ampl_constraint_index) {
		if (PRINT_) std::cout << "On constraint ampl_constraint_index=" << ampl_constraint_index << ": ";
		bool equality_con = LUrhs[2*ampl_constraint_index]==LUrhs[2*ampl_constraint_index+1] && LUrhs[2*ampl_constraint_index]!=INFINITY;
	
		if (equality_con) {
			if (PRINT_) std::cout << "equality";
			if (PRINT_) std::cout << "(push eq)";
			equality_constraints_.push_back(ampl_constraint_index);
		} else if (LUrhs[2*ampl_constraint_index] == -INFINITY && LUrhs[2*ampl_constraint_index+1] > -INFINITY) {
			if (PRINT_) std::cout << "upper bound";
			if (PRINT_) std::cout << "(push ineq)";
			inequality_constraints_upper_.push_back(ampl_constraint_index);
		} else if (LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] < INFINITY) {
			if (PRINT_) std::cout << "double-sided";
			if (PRINT_) std::cout << "(push ineq)";
			if (PRINT_) std::cout << "(push -ineq)";
			inequality_constraints_upper_.push_back(ampl_constraint_index);
			inequality_constraints_lower_.push_back(ampl_constraint_index);
		} else if (LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] == INFINITY) {
			if (PRINT_) std::cout << "lower bound only";
			if (PRINT_) std::cout << "(push -ineq)";
			inequality_constraints_lower_.push_back(ampl_constraint_index);
		} else {
			std::cout << "unsupported";
		}
		if (PRINT_) std::cout << std::endl;
	}
	
	
	
	std::vector<double> x(this->num_primal());
	// double *x = new double[num_primal()];
	for (size_t primal_index=0; primal_index < this->num_primal(); ++primal_index)
		x[primal_index] = X0_[primal_index];
	
	bool should_abort = false;
	for (size_t ampl_variable_index =0 ; ampl_variable_index < n_var; ++ampl_variable_index) {
		if (PRINT_) 
			std::cout << "ampl_variable_index=" << ampl_variable_index << ": l=" << LUv[2*ampl_variable_index] << " <= x_ampl_variable_index= " << x[ampl_variable_index] << " <= u=" << LUv[2*ampl_variable_index+1] << std::endl;
		bool lower_finite = (LUv[2*ampl_variable_index] != -INFINITY);
		bool upper_finite = (LUv[2*ampl_variable_index+1] != INFINITY);
		bool both_finite = lower_finite && upper_finite;
		bool equal = (LUv[2*ampl_variable_index] == LUv[2*ampl_variable_index+1]);
		if (both_finite && equal) {
			variable_equality_.push_back(ampl_variable_index);
		} else if (both_finite) {
			variable_bound_lower_.push_back(ampl_variable_index);
			variable_bound_upper_.push_back(ampl_variable_index);
		} else if (lower_finite) {
			variable_bound_lower_.push_back(ampl_variable_index);
		} else if (upper_finite) {
			variable_bound_upper_.push_back(ampl_variable_index);
		}			
	}
	
	std::vector<double> con(n_con);
    num_constraint_eval_++;
	conval(&x[0], &con[0], nerror_);
	for (size_t ampl_constraint_index =0 ; ampl_constraint_index < n_con; ++ampl_constraint_index) {
		if (PRINT_) std::cout << "ampl_constraint_index=" << ampl_constraint_index << ": l=" << LUrhs[2*ampl_constraint_index] << " <= c(x_k)= " << con[ampl_constraint_index] << " <= u=" << LUrhs[2*ampl_constraint_index+1] << std::endl;
	}
	
	this->num_dual_eq_ = equality_constraints_.size();
	this->num_dual_ieq_ = inequality_constraints_lower_.size() + inequality_constraints_upper_.size() + variable_bound_lower_.size() + variable_bound_upper_.size();
	asl_ = asl;
}

AmplNlp::~AmplNlp() {
	ASL *asl = asl_;
	
	delete nerror_;
	ASL_free(&asl);
	fclose(nl_);
}

iSQOIterate AmplNlp::initial(double penalty_parameter) {
    iSQOIterate initial(this->num_primal(), this->num_dual_eq(), this->num_dual_ieq(), penalty_parameter);
    initial.assign_primal(X0_);

	std::vector<double> eq_values = constraints_equality(initial);
    std::vector<double> dual_eq_values(this->num_dual_eq());
	for (size_t dual_eq_index=0; dual_eq_index < this->num_dual_eq(); ++dual_eq_index) {
		if (eq_values[dual_eq_index] < 0) {
			dual_eq_values[dual_eq_index] = -1.0;
		} else if (eq_values[dual_eq_index] > 0) {
			dual_eq_values[dual_eq_index] = +1.0;
		} else {
			dual_eq_values[dual_eq_index] = 0.0;
		}
	}
    initial.assign_dual_eq(dual_eq_values);

	std::vector<double> ieq_values = constraints_inequality(initial);
    std::vector<double> dual_ieq_values(this->num_dual_ieq());
	for (size_t dual_ieq_index=0; dual_ieq_index < this->num_dual_ieq(); ++dual_ieq_index) {
		if (ieq_values[dual_ieq_index] < 0) {
			dual_ieq_values[dual_ieq_index] = 0.0;
		} else if (ieq_values[dual_ieq_index] > 0) {
			dual_ieq_values[dual_ieq_index] =+1.0;
		} else {
			dual_ieq_values[dual_ieq_index] = 0.0;
		}
	}
    initial.assign_dual_ieq(dual_ieq_values);
    
	return initial;
}
// zeroth order NLP quantities:
double AmplNlp::objective(const iSQOIterate &iterate) {
	ASL *asl = asl_;
    // TODO look into doing const casts here (& in all other AmplNlp functions...), so that we don't have to do any copies!
    std::vector<double> x(this->num_primal());
    x.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());

    num_objective_eval_++;
    real obj = objval(0, &x[0], nerror_);
	if (PRINT_) std::cout << "objective value(nerror = " << *nerror_ << "): " << obj << std::endl;
	return obj;
}

std::vector<double> AmplNlp::constraints_equality(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(this->num_primal());
    x.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
    
	std::vector<double> con(n_con);
    num_constraint_eval_++;
	conval(&x[0], &con[0], nerror_);
	if (PRINT_) std::cout << "equality_constraints:" << std::endl;
	std::vector<double> equality_constraint_evaluation(equality_constraints_.size());
	for (size_t isqo_eq_constraint_index=0; isqo_eq_constraint_index < equality_constraints_.size(); ++isqo_eq_constraint_index) {
		equality_constraint_evaluation[isqo_eq_constraint_index] = con[equality_constraints_[isqo_eq_constraint_index]] - LUrhs[2*equality_constraints_[isqo_eq_constraint_index]];
	}
	if (PRINT_) std::cout << "eq eval: " << equality_constraint_evaluation;
	return equality_constraint_evaluation;
}

std::vector<double> AmplNlp::constraints_inequality(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(this->num_primal());
    x.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
    
    if (PRINT_) std::cout << "my value of x in AmplNlp::constraints_inequality: " << *(iterate.get_primal_values()) << std::endl;
    if (PRINT_) std::cout << "my value of x in AmplNlp::constraints_inequality: " << x << std::endl;
    std::vector<double> con(n_con);
	
    num_constraint_eval_++;
	conval(&x[0], &con[0], nerror_);
	if (PRINT_) std::cout << "inequality_constraints:" << std::endl;
	
	std::vector<double> inequality_constraint_evaluation(this->num_dual_ieq());
	
	int isqo_ineq_constraint_index = 0;
	for (size_t isqo_ineq_lower_constraint_index=0; isqo_ineq_lower_constraint_index < inequality_constraints_lower_.size(); ++isqo_ineq_lower_constraint_index) {
		if (PRINT_) std::cout << "isqo ineq " << isqo_ineq_constraint_index << "; " << inequality_constraints_lower_[isqo_ineq_lower_constraint_index] << std::endl;
		size_t ampl_constraint_index = inequality_constraints_lower_[isqo_ineq_lower_constraint_index];
		size_t ampl_constraint_value_lower = 2*ampl_constraint_index;
		inequality_constraint_evaluation[isqo_ineq_constraint_index] = -con[ampl_constraint_index] + LUrhs[ampl_constraint_value_lower];
		++isqo_ineq_constraint_index;
	}
	for (size_t isqo_ineq_upper_constraint_index=0; isqo_ineq_upper_constraint_index < inequality_constraints_upper_.size(); ++isqo_ineq_upper_constraint_index) {
		
		size_t ampl_constraint_index = inequality_constraints_upper_[isqo_ineq_upper_constraint_index];
		size_t	ampl_constraint_value_upper = 2*ampl_constraint_index+1;
		if (PRINT_) std::cout << "isqo ineq " << isqo_ineq_upper_constraint_index << " is AMPL constraint " << ampl_constraint_index << std::endl;
		inequality_constraint_evaluation[isqo_ineq_constraint_index] = con[ampl_constraint_index] - LUrhs[ampl_constraint_value_upper];
		++isqo_ineq_constraint_index;
	}
	for (size_t isqo_ineq_lower_variable_index=0; isqo_ineq_lower_variable_index < variable_bound_lower_.size(); ++isqo_ineq_lower_variable_index) {
		inequality_constraint_evaluation[isqo_ineq_constraint_index] = LUv[2*variable_bound_lower_[isqo_ineq_lower_variable_index]] - x[variable_bound_lower_[isqo_ineq_lower_variable_index]];
		// std::cout 
		// 	<< "variable lower bound " << isqo_ineq_lower_variable_index 
		// 	<< "; ampl bound: " <<  variable_bound_lower_[isqo_ineq_lower_variable_index]
		// 	<< "; with value: " << LUv[2*variable_bound_lower_[isqo_ineq_lower_variable_index]]
		// 	<< "; @x_k: " << LUv[2*variable_bound_lower_[isqo_ineq_lower_variable_index]] - x[variable_bound_lower_[isqo_ineq_lower_variable_index]]
		// 	<< std::endl;
		++isqo_ineq_constraint_index;
	}
	for (size_t isqo_ineq_upper_variable_index=0; isqo_ineq_upper_variable_index < variable_bound_upper_.size(); ++isqo_ineq_upper_variable_index) {
		// std::cout << "variable upper bound " << isqo_ineq_upper_variable_index << "" << std::endl;
        // POTENTIALLY BUGGY LINE RIGHT HERE:
        // inequality_constraint_evaluation[isqo_ineq_constraint_index] = -LUv[2*variable_bound_lower_[isqo_ineq_upper_variable_index]+1] + x[variable_bound_lower_[isqo_ineq_upper_variable_index]];
        inequality_constraint_evaluation[isqo_ineq_constraint_index] = -LUv[2*variable_bound_upper_[isqo_ineq_upper_variable_index]+1] + x[variable_bound_upper_[isqo_ineq_upper_variable_index]];
		++isqo_ineq_constraint_index;
	}
	assert(isqo_ineq_constraint_index == this->num_dual_ieq());
	
	
	if (PRINT_) std::cout <<"ineq eval" << inequality_constraint_evaluation;
	return inequality_constraint_evaluation;
}

// first order NLP quantities:
std::vector<double> AmplNlp::objective_gradient(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> return_gradient(this->num_primal());
	
	std::vector<double> x(this->num_primal());
    x.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
    if (PRINT_) std::cout << "evaluating the gradient with this iterate: " << iterate << std::endl;
    if (PRINT_) std::cout << "evaluating the gradient with this x: " << x << std::endl;
	num_objective_gradient_eval_++;
	objgrd(0, &x[0], &return_gradient[0], nerror_);
	if (PRINT_) std::cout << "objective gradient(nerror = " << *nerror_ << "): " << return_gradient[0] << std::endl;
	if (PRINT_) std::cout << "objective gradient(nerror = " << *nerror_ << "): " << return_gradient[1] << std::endl;
    if (PRINT_) std::cout << "AmplNlp::objective_gradient return_gradient : " << return_gradient << std::endl;
	
	return return_gradient;
}

std::shared_ptr<matrix_base_class> DenseAmplNlp::constraints_equality_jacobian(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(this->num_primal());
    x.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
	
    // std::shared_ptr<matrix_base_class> equality_constraint_jacobian(new dense_matrix(equality_constraints_.size(), n_var));
    std::shared_ptr<dense_matrix> equality_constraint_jacobian = std::shared_ptr<dense_matrix>(new dense_matrix(equality_constraints_.size(), this->num_primal()));
	if (PRINT_) std::cout << "equal: " << equality_constraints_ << std::endl;
	std::vector<double> G(num_primal());
    num_constraint_jacobian_eval_++;
	for (size_t isqo_eq_constraint_index=0; isqo_eq_constraint_index < equality_constraints_.size(); ++isqo_eq_constraint_index) {
		congrd(equality_constraints_[isqo_eq_constraint_index], &x[0], &G[0], nerror_);
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			equality_constraint_jacobian->set(isqo_eq_constraint_index, primal_index, G[primal_index]);
		}
		
	}
	if (PRINT_) std::cout << "equality jacobian: [" << std::endl;
	if (PRINT_) std::cout << equality_constraint_jacobian;
	if (PRINT_) std::cout << "]" << std::endl;

	return std::shared_ptr<matrix_base_class>(equality_constraint_jacobian);
}

std::shared_ptr<sparse_matrix> SparseAmplNlp::submatrix(std::shared_ptr<sparse_matrix> matrix, std::vector<size_t> rows_to_keep, size_t num_result_nnz, double scale) {
    
    // create a keep_rows vector so that we know which rows to keep(avoid O(n) searches through constraints index.)
    std::vector<int> keep_rows(matrix->num_rows(), -1);
    if (PRINT_) std::cout << "keep_rows: " << matrix->num_columns() << "; " << keep_rows << std::endl;
    int current_result_row_tally=0;
    for (size_t keep_rows_index=0; keep_rows_index < rows_to_keep.size(); ++keep_rows_index) {
        keep_rows[rows_to_keep[keep_rows_index]] = current_result_row_tally++;
        
    }
    if (PRINT_) std::cout << "keep_rows: " << matrix->num_columns() << "; " << keep_rows << std::endl;
    
    std::shared_ptr<sparse_matrix> result = std::shared_ptr<sparse_matrix>(new sparse_matrix(rows_to_keep.size(), matrix->num_columns(), num_result_nnz));
    if (PRINT_) std::cout << "result: " << result << std::endl;
    
    std::vector<int> result_col_counts(matrix->num_columns()+1);
    int current_result_nonzero = 0;
    for (size_t col_index=0; col_index < matrix->num_columns(); ++col_index) {
        for (size_t master_nonzero_index=matrix->get_col_start(col_index); master_nonzero_index < matrix->get_col_start(col_index+1); ++master_nonzero_index) {
            if (keep_rows[matrix->get_row_index(master_nonzero_index)] >= 0) {
                result->set_value(current_result_nonzero, scale * matrix->get_value(master_nonzero_index));
                result->set_row_index(current_result_nonzero, keep_rows[matrix->get_row_index(master_nonzero_index)]);
                ++result_col_counts[col_index];
                ++current_result_nonzero;
            }
        }
        result->set_col_start(col_index+1, current_result_nonzero);
    }
    if (PRINT_) std::cout << "result: " << result << std::endl;
    
    // assert(false);
    return result;
}

void SparseAmplNlp::jacobian_update(const iSQOIterate &iterate) {
	ASL *asl = asl_;
    PRINT_ = false;
	std::vector<double> x(num_primal());
    x.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
		
    asl->i.congrd_mode = 0;
    double J[nzc];
    jacval(&x[0], &J[0], nerror_);
    num_constraint_jacobian_eval_++;
    if (PRINT_) std::cout << "J: " << std::endl;
    for (size_t current_nz=0; current_nz < nzc; ++current_nz) {
        if (PRINT_) std::cout << "current_nz=" << current_nz << ", val=" << J[current_nz] << std::endl;
    }
    
    cgrad *cg;
    size_t my_count_nz =0;
    
    if (PRINT_) std::cout << "equ: " << equality_constraints_ << std::endl;
    if (PRINT_) std::cout << "ineq lower: " << inequality_constraints_lower_ << std::endl;
    if (PRINT_) std::cout << "ineq upper: " << inequality_constraints_upper_ << std::endl;
        
    if (PRINT_) std::cout << "FULL TABLE: " << std::endl;
    // the number of nonzero entries in each type:
    size_t num_eq_nnz=0, num_ieq_lower_nnz=0, num_ieq_upper_nnz=0;
    std::shared_ptr<sparse_matrix> ampl_full_jacobian = std::shared_ptr<sparse_matrix>(new sparse_matrix(n_con, n_var, nzc));
    std::vector<int> col_counts(n_var+1);
    for (size_t ampl_constraint_index=0; ampl_constraint_index<n_con; ++ampl_constraint_index) {
        for (cg = Cgrad[ampl_constraint_index]; cg; cg = cg->next){
            char constraint_type;
            if ( LUrhs[2*ampl_constraint_index] == LUrhs[2*ampl_constraint_index+1] ) {
                constraint_type = 'E';
                ++num_eq_nnz;
            } else if ( LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] < INFINITY ) { 
                constraint_type = 'B';
                ++num_ieq_lower_nnz;
                ++num_ieq_upper_nnz;
            } else if ( LUrhs[2*ampl_constraint_index] == -INFINITY && LUrhs[2*ampl_constraint_index+1] < INFINITY ) { 
                constraint_type = 'U';
                ++num_ieq_upper_nnz;
            } else if ( LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] == INFINITY ) { 
                constraint_type = 'L';
                ++num_ieq_lower_nnz;
            } else {
                assert(false);
            }
            
            if (PRINT_) printf("index: %3d: conindex: %3d[%9.2e <= c_i <= %9.2e ==> %c]: varno: %3d, coef: %9.2e, goff: %3d, val: %9.2e\n", 
                    my_count_nz, ampl_constraint_index, LUrhs[2*ampl_constraint_index], LUrhs[2*ampl_constraint_index+1], constraint_type,
                    cg->varno, cg->coef, cg->goff, J[cg->goff]);

            ampl_full_jacobian->set_value(cg->goff, J[cg->goff]);
            ampl_full_jacobian->set_row_index(cg->goff, ampl_constraint_index);
            ++col_counts[cg->varno];
            ++my_count_nz;
        }
    }
    if (PRINT_) std::cout << "col counts: " << col_counts << std::endl;
    std::vector<int> col_starts(n_var+1);
    int accumulate_num_nnz=0;
    for (size_t col_index=0; col_index < n_var; ++col_index) {
        accumulate_num_nnz += col_counts[col_index];
        // col_starts[col_index+1] = accumulate_num_nnz;
        ampl_full_jacobian->set_col_start(col_index+1, accumulate_num_nnz);
    }
    
    
    if (PRINT_) std::cout << "ampl_full_jacobian: " << ampl_full_jacobian << std::endl;    
    if (PRINT_) std::cout << "=================================================== (have full AMPL jacobian now.)" << std::endl;
    
    eq_jacobian_ = submatrix(ampl_full_jacobian, equality_constraints_, num_eq_nnz, 1.0);
    std::shared_ptr<sparse_matrix> ieq_lower_jacobian = submatrix(ampl_full_jacobian, inequality_constraints_lower_, num_ieq_lower_nnz, -1.0);
    std::shared_ptr<sparse_matrix> ieq_upper_jacobian = submatrix(ampl_full_jacobian, inequality_constraints_upper_, num_ieq_upper_nnz, +1.0);
    
    if (PRINT_) std::cout << "eq_jacobian_: " << eq_jacobian_ << std::endl;
    if (PRINT_) std::cout << "ieq_lower_jacobian: " << ieq_lower_jacobian << std::endl;
    if (PRINT_) std::cout << "ieq_upper_jacobian: " << ieq_upper_jacobian << std::endl;
    
    if (PRINT_) std::cout << "constraint inequality ------------" << std::endl;
    std::shared_ptr<sparse_matrix> constraints_ieq = ieq_lower_jacobian->vertical(ieq_upper_jacobian);
    if (PRINT_) std::cout << constraints_ieq << std::endl;

    std::shared_ptr<sparse_matrix> lower_var(new sparse_matrix(n_var, variable_bound_lower_, -1.0));
    std::shared_ptr<sparse_matrix> upper_var(new sparse_matrix(n_var, variable_bound_upper_, +1.0));
    if (PRINT_) std::cout << "lower_var: " << lower_var << std::endl;
    if (PRINT_) std::cout << "upper_var: " << upper_var << std::endl;
    std::shared_ptr<sparse_matrix> variable_bounds_jacobian = lower_var->vertical(upper_var);
    if (PRINT_) std::cout << "variable_bounds_jacobian: " << variable_bounds_jacobian << std::endl;
    
    if (PRINT_) std::cout << "=========\n\n\n" << std::endl;
    
    
    ieq_jacobian_ = constraints_ieq->vertical(variable_bounds_jacobian);
    if (PRINT_) std::cout << "full inequality ------------" << std::endl;
    if (PRINT_) std::cout << ieq_jacobian_ << std::endl;
    
    return;
}

std::shared_ptr<matrix_base_class> SparseAmplNlp::constraints_equality_jacobian(const iSQOIterate &iterate) {
    std::map<int, std::shared_ptr<matrix_base_class> >::iterator it = prior_eq_jacobians_.find(iterate.get_serial());
    if (it != prior_eq_jacobians_.end()) return it->second;
    
    jacobian_update(iterate);
    prior_eq_jacobians_[iterate.get_serial()] = eq_jacobian_;
	return eq_jacobian_;
}

std::shared_ptr<matrix_base_class> SparseAmplNlp::constraints_inequality_jacobian(const iSQOIterate &iterate) {
    std::map<int, std::shared_ptr<matrix_base_class> >::iterator it = prior_ieq_jacobians_.find(iterate.get_serial());
    if (it != prior_ieq_jacobians_.end()) return it->second;
    
    jacobian_update(iterate);
    prior_ieq_jacobians_[iterate.get_serial()] = ieq_jacobian_;
	return ieq_jacobian_;
}

std::shared_ptr<matrix_base_class> DenseAmplNlp::constraints_inequality_jacobian(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(num_primal());
    x.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
    
    // std::shared_ptr<matrix_base_class> inequality_constraint_jacobian(new dense_matrix(num_dual_ieq(), n_var));
    dense_matrix *inequality_constraint_jacobian= new dense_matrix(num_dual_ieq(), num_primal());
	size_t isqo_ineq_constraint_index=0;
	if (PRINT_) std::cout << "lower: " << inequality_constraints_lower_ << std::endl;
	if (PRINT_) std::cout << "upper: " << inequality_constraints_upper_ << std::endl;
	std::vector<double> G(num_primal());
	for (size_t isqo_ineq_lower_constraint_index=0; isqo_ineq_lower_constraint_index < inequality_constraints_lower_.size(); ++isqo_ineq_lower_constraint_index) {
		congrd(inequality_constraints_lower_[isqo_ineq_lower_constraint_index], &x[0], &G[0], nerror_);
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			inequality_constraint_jacobian->set(isqo_ineq_constraint_index, primal_index, -G[primal_index]);
		}
		++isqo_ineq_constraint_index;
	}
	for (size_t isqo_ineq_upper_constraint_index=0; isqo_ineq_upper_constraint_index < inequality_constraints_upper_.size(); ++isqo_ineq_upper_constraint_index) {
		congrd(inequality_constraints_upper_[isqo_ineq_upper_constraint_index], &x[0], &G[0], nerror_);
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			inequality_constraint_jacobian->set(isqo_ineq_constraint_index, primal_index, G[primal_index]);
		}
		++isqo_ineq_constraint_index;
	}
    // set up Jacobian for variable bounds.
    // LOWER VARIABLE BOUNDS:
	for (size_t isqo_ineq_lower_variable_index=0; isqo_ineq_lower_variable_index < variable_bound_lower_.size(); ++isqo_ineq_lower_variable_index) {
		inequality_constraint_jacobian->set(isqo_ineq_constraint_index, variable_bound_lower_[isqo_ineq_lower_variable_index], -1.0);
		++isqo_ineq_constraint_index;
	}
    // UPPER VARIABLE BOUNDS:
	for (size_t isqo_ineq_upper_variable_index=0; isqo_ineq_upper_variable_index < variable_bound_upper_.size(); ++isqo_ineq_upper_variable_index) {
		inequality_constraint_jacobian->set(isqo_ineq_constraint_index, variable_bound_upper_[isqo_ineq_upper_variable_index], +1.0);
		++isqo_ineq_constraint_index;
	}
    
	if (PRINT_) 
	{
		std::cout << "======" << std::endl;
		std::cout << inequality_constraint_jacobian << std::endl;
		std::cout << "======" << std::endl;
	}
			
	return std::shared_ptr<matrix_base_class>(inequality_constraint_jacobian);
}

std::vector<double> AmplNlp::mux_multipliers(const iSQOIterate &iterate) {
    ASL *asl = asl_;
	std::vector<double> Y(n_con);
	
	// funnel dual_{,i}eq_values_ into Y to pass to AMPL.
    std::vector<double> dual_eq_values(this->num_dual_eq());
    dual_eq_values.assign(iterate.get_dual_eq_values()->begin(), iterate.get_dual_eq_values()->end());
	for (size_t dual_eq_index=0; dual_eq_index < equality_constraints_.size(); ++dual_eq_index) {
		Y[equality_constraints_[dual_eq_index]] = dual_eq_values[dual_eq_index];// / iterate.penalty_parameter_;
	}
	
    std::vector<double> dual_ieq_values(this->num_dual_ieq());
    dual_ieq_values.assign(iterate.get_dual_ieq_values()->begin(), iterate.get_dual_ieq_values()->end());
	size_t dual_ineq_value_index = 0;
	for (size_t dual_ineq_lower_index=0; dual_ineq_lower_index < inequality_constraints_lower_.size(); ++dual_ineq_lower_index) {
		Y[inequality_constraints_lower_[dual_ineq_lower_index]] = -dual_ieq_values[dual_ineq_value_index];
		dual_ineq_value_index++;
	}
	for (size_t dual_ineq_upper_index=0; dual_ineq_upper_index < inequality_constraints_upper_.size(); ++dual_ineq_upper_index) {
		Y[inequality_constraints_upper_[dual_ineq_upper_index]] = dual_ieq_values[dual_ineq_value_index];
		dual_ineq_value_index++;
	}
    return Y;
}
// second order NLP quantities:
std::shared_ptr<matrix_base_class> DenseAmplNlp::lagrangian_hessian(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> H(num_primal() * num_primal());
	std::vector<double> OW(1);
	OW[0] = iterate.get_penalty_parameter();
	std::vector<double> Y = mux_multipliers(iterate);
	
	// In this call:
	//	 - H : OUTPUT : vector holding entries of full hessian.
	//	 - n_var : INPUT : # of variables or something about the stride, maybe.
	//	 - 0 : INPUT : index of the desired objective.
	//	 - OW : INPUT : multipliers for objective function ("objective weights")
	//	 - Y : INPUT : lagrange multipliers for constraints.
	fullhes(&H[0], n_var, 0, &OW[0], &Y[0]);
    num_lagrangian_eval_++;
    // the code below is highly suspect.
    // I think it was necesary in the MATLAB implementation, where we couldn't set
    // the objective multiplier, but since we can control that with OW[0] now, 
    // there's no reason not to.
	// if (iterate.penalty_parameter_ == 0) {
//         std::vector<double> H_f_only(num_primal()*num_primal());
//         std::vector<double> OW_f_only(1);
//         OW_f_only[0] = 1.0;
//         std::vector<double> Y_f_only(n_con);
//         for (size_t dual_index=0; dual_index < n_con; ++dual_index) Y[dual_index] = 0.0;
//         fullhes(&H_f_only[0], n_var, 0, &OW_f_only[0], &Y_f_only[0]);
//         
//         for (size_t row_index=0; row_index < n_var; ++row_index)
//             for (size_t col_index=0; col_index < n_var; ++col_index) 
//                 H[row_index*n_var + col_index] -= H_f_only[row_index*n_var + col_index];
//     }
	
    dense_matrix *return_hessian = new dense_matrix(num_primal(), num_primal());
    // std::shared_ptr<matrix_base_class> return_hessian(new dense_matrix(n_var,n_var));
	for (size_t row_index=0; row_index < n_var; ++row_index) {
		for (size_t col_index=0; col_index < n_var; ++col_index) {
			return_hessian->set(row_index, col_index, H[row_index*n_var + col_index]);
		}
	}
	return std::shared_ptr<matrix_base_class>(return_hessian);
}

// void AmplNlp::hessian_update(const iSQOIterate &iterate) {
    // 
// }

std::shared_ptr<matrix_base_class> SparseAmplNlp::lagrangian_hessian(const iSQOIterate &iterate) {

    if (PRINT_) std::cout << "trying to create a sparse hessian using iterate serial #: " << iterate.get_serial() << std::endl;
    
    if (PRINT_) std::cout << "testing for a cache hit." << std::endl;
    std::map<int, std::shared_ptr<matrix_base_class> >::iterator it;
    size_t index = 0;
    for (std::map<int, std::shared_ptr<matrix_base_class>>::iterator it = prior_hessians_.begin(); it != prior_hessians_.end(); it++) {
        if (PRINT_) std::cout << "testing index=" << index << "; " << it->first << ", "; // << *(it->second)
        if (it->first == iterate.get_serial()) {
            if (PRINT_) std::cout << "HIT HIT HIT HIT" << std::endl;
            return it->second;
        } else {
            if (PRINT_) std::cout << "MISS MISS MISS" << std::endl;
        }
        index++;
    }
    
    
    PRINT_ = false;
    ASL *asl = asl_;
	std::vector<double> H(num_primal() * num_primal());
	std::vector<double> OW(1);
	OW[0] = iterate.get_penalty_parameter();
	std::vector<double> Y = mux_multipliers(iterate);
	
	// In this call:
	//	 - H : OUTPUT : vector holding entries of full hessian.
	//	 - n_var : INPUT : # of variables or something about the stride, maybe.
	//	 - 0 : INPUT : index of the desired objective.
	//	 - OW : INPUT : multipliers for objective function ("objective weights")
	//	 - Y : INPUT : lagrange multipliers for constraints.
	// fullhes(&H[0], n_var, 0, &OW[0], &Y[0]);
//     if (PRINT_) {
//         for (size_t i=0; i<n_var; ++i) {
//             for (size_t j=0; j<n_var; ++j) {
//                 std::cout << "H[n_var*i=" << i << ", j=" << j << "]: " << H[n_var*i + j] << std::endl;
//             }
//         }
//     }
    
    
    // fint sphsetup(int nobj, int ow, int y, int uptri)
    int hessian_nnz = sphsetup(-1, 1, 1, 0);
    
    if (PRINT_) {
        std::cout << "nonzeros in hessian: " << hessian_nnz << std::endl;
    
        for (size_t current_nonzero=0; current_nonzero<hessian_nnz; ++current_nonzero) {
            std::cout << "sputinfo->hrownos[current_nonzero = " << current_nonzero << "]: " << sputinfo->hrownos[current_nonzero] << std::endl;
        }
        for (size_t current_column=0; current_column<n_var+1; ++current_column) {
            std::cout << "sputinfo->hcolstarts[current_column = " << current_column<< "]: " << sputinfo->hcolstarts[current_column] << std::endl;
        }
    }
    
    std::vector<double> sparse_hessian_values(hessian_nnz);
    // void sphes(real *H, int nobj, real *OW, real *Y)
    // sphes arguments:
    //  - 0: *H: pointer to array to fill nonzero values into.
    //  - 1: nobj: which objective function to use (-1 specifies 'all', counter to all the documentation & dense syntax.)
    //  - 2: *OW: pointer to array of 'objective weights'
    //  - 3: *Y: pointer to array of lagrange multipliers.
    sphes(&sparse_hessian_values[0], -1, &OW[0], &Y[0]); // &OW[0], &Y[0]
    
    num_lagrangian_eval_++;
    if (PRINT_) 
        for (size_t sparse_hessian_index=0; sparse_hessian_index<hessian_nnz; ++sparse_hessian_index) {
            std::cout << "sparse_hessian_values[sparse_hessian_index=" << sparse_hessian_index << "]: " << sparse_hessian_values[sparse_hessian_index] << std::endl;
        }
    
    if (PRINT_) std::cout << "===========" << std::endl;
    // std::shared_ptr<matrix_base_class> sparse_hessian_(new sparse_matrix(n_var, n_var, hessian_nnz));
    sparse_matrix *sparse_hessian_ = new sparse_matrix(n_var, n_var, hessian_nnz);
    std::vector<int> sparse_hessian_col_counts(n_var);
    std::vector<int> sparse_hessian_accumulation_counter(n_var);
    std::vector<int> sparse_hessian_accumulation(hessian_nnz);
    std::vector<int> ampl_column_ind(hessian_nnz);
    for (size_t current_column=0; current_column < n_var; ++current_column) {
        if (PRINT_) std::cout << "col: " << current_column << ": nonzeros " << sputinfo->hcolstarts[current_column] << " to " << sputinfo->hcolstarts[current_column+1] << std::endl;
        for (size_t current_nonzero=sputinfo->hcolstarts[current_column]; current_nonzero < sputinfo->hcolstarts[current_column+1]; ++current_nonzero) {
            int ampl_row = sputinfo->hrownos[current_nonzero];
            int ampl_col = current_column;
            int qpOASES_row = ampl_col;
            int qpOASES_col = ampl_row;
            
            if (PRINT_) std::cout << " +-> current_nonzero: " << current_nonzero << "; row: " << ampl_row << std::endl;
            // std::cout << " ---- new entry: " << 
            ++sparse_hessian_col_counts[qpOASES_col];
            sparse_hessian_accumulation[current_nonzero] = sparse_hessian_accumulation_counter[ampl_row]; // technically unnecessary...
            ampl_column_ind[current_nonzero] = current_column;
            ++sparse_hessian_accumulation_counter[ampl_row];
        }
    }
    
    if (PRINT_) std::cout << "sparse_hessian_col_counts  : " << sparse_hessian_col_counts << std::endl;
    if (PRINT_) std::cout << "sparse_hessian_accumulation_counter: " << sparse_hessian_accumulation_counter << std::endl;
    if (PRINT_) std::cout << "sparse_hessian_accumulation: " << sparse_hessian_accumulation << std::endl;
    int current_nonzero_count=0;
    for (size_t current_column=0; current_column < n_var; ++current_column) {
        current_nonzero_count += sparse_hessian_col_counts[current_column];
        sparse_hessian_->set_col_start(current_column+1, current_nonzero_count);
    }
    // if (PRINT_) std::cout << "col starts: " << sparse_hessian_->col_starts_ << std::endl;
    
    if (PRINT_) std::cout << "===========" << std::endl;
    std::vector<int> ampl_to_qp_indices(hessian_nnz);
    for (size_t current_column=0; current_column < n_var; ++current_column) {
        if (PRINT_) std::cout << "col: " << current_column << ": nonzeros " << sputinfo->hcolstarts[current_column] << " to " << sputinfo->hcolstarts[current_column+1] << std::endl;
        for (size_t current_nonzero=sputinfo->hcolstarts[current_column]; current_nonzero < sputinfo->hcolstarts[current_column+1]; ++current_nonzero) {
            int ampl_row = sputinfo->hrownos[current_nonzero];
            int ampl_col = current_column;
            int qpOASES_row = ampl_col;
            int qpOASES_col = ampl_row;
            
            ampl_to_qp_indices[current_nonzero] = sparse_hessian_->get_col_start(ampl_row) + sparse_hessian_accumulation[current_nonzero];
            
            if (PRINT_) std::cout << " +-> current_nonzero: " << current_nonzero << "; "
                        << "row: " << ampl_row << "; "
                        << "qp_col_st[ampl_row]: " << sparse_hessian_->get_col_start(ampl_row) << "; "
                        << "sparse_hessian_accumulation[current_nonzero]: " << sparse_hessian_accumulation[current_nonzero]
                        // << std::endl
                            ;
            if (PRINT_) std::cout // << " +--> "
                << "; ampl_to_qp_indices: " << ampl_to_qp_indices[current_nonzero]
                << "; ampl_cols[<-]: " << ampl_column_ind[ampl_to_qp_indices[current_nonzero]]
                << std::endl;
            
            sparse_hessian_->set_value(current_nonzero, sparse_hessian_values[ampl_to_qp_indices[current_nonzero]]);
            sparse_hessian_->set_row_index(current_nonzero, ampl_column_ind[ampl_to_qp_indices[current_nonzero]]);
        }
    }
    if (PRINT_) std::cout << "ampl_to_qp_indices: " << ampl_to_qp_indices << std::endl;
    
    if (PRINT_) std::cout << "sparse_hessian_: " << sparse_hessian_ << std::endl;
    std::shared_ptr<matrix_base_class> retval(sparse_hessian_);
    prior_hessians_[iterate.get_serial()] = retval;
    return retval;
}


std::size_t AmplNlp::num_lower_ieq() {
	return inequality_constraints_lower_.size();
}
std::size_t AmplNlp::num_upper_ieq() {
	return this->inequality_constraints_upper_.size();
}

