#ifndef ISQO_CONFIG_CC_3N2PPVP3
#define ISQO_CONFIG_CC_3N2PPVP3

#include <fstream>
#include <iostream>
#include <map>
// initial_penalty_parameter: [numerical]
// hessian_shift_print_level: none | full
// qp_print_level: none | min | medium | max | table
// termination_tolerance: [numerical]
// 
enum hessian_shift_print_level_t {HSPL_NONE, HSPL_FULL};
enum qp_print_level_t {QPL_NONE, QPL_LOW, QPL_MEDIUM, QPL_HIGH, QPL_TABLE};

class iSQOControlPanel {
public:    
    // TODO should be a singleton??
    iSQOControlPanel(std::string config_filename);
    
    template <typename T, typename S>
    void update_junk(std::istream &is, std::map<std::string, T> update_map, S &value_to_update);
    void update_settings();
    
    void print();
    hessian_shift_print_level_t get_hessian_shift_print_level();
    double get_initial_penalty_parameter();
    qp_print_level_t get_qp_print_level();
    double get_termination_tolerance();
    
private:
    hessian_shift_print_level_t hessian_shift_print_level;
    double initial_penalty_parameter;
    qp_print_level_t qp_print_level;
    double termination_tolerance;
    
    std::map<std::string, hessian_shift_print_level_t> map_hessian_shift_print_levels;
    std::map<std::string, qp_print_level_t> map_qp_print_levels;
    
    std::ifstream config_file_;
};


#endif /* end of include guard: ISQO_CONFIG_CC_3N2PPVP3 */
