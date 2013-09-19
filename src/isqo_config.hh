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
    
    typedef void (*UpdateFunction)(void);
    
    // TODO should be a singleton??
    iSQOControlPanel() {
        hessian_shift_print_level = HSPL_NONE;
        initial_penalty_parameter = 1e-1;
        qp_print_level=QPL_NONE;
        termination_tolerance=1e-6;
        
        // setup the maps from config strings -> *_t
        map_hessian_shift_print_levels["none"] = HSPL_NONE;
        map_hessian_shift_print_levels["full"] = HSPL_FULL;
        map_qp_print_levels["none"] = QPL_NONE;
        map_qp_print_levels["low"] = QPL_LOW;
        map_qp_print_levels["medium"] = QPL_MEDIUM;
        map_qp_print_levels["high"] = QPL_HIGH;
        map_qp_print_levels["table"] = QPL_TABLE;
    } 
    
    template <typename T, typename S>
    void update_junk(std::istream &is, std::map<std::string, T> update_map, S &value_to_update) {
        std::string scratch;
        is >> scratch;
        if (update_map.find(scratch) != update_map.end()) {
            value_to_update = update_map[scratch];
        } else {
            std::cerr << "unknown " << "??? option specified: [" << scratch << "]."<< std::endl;
        }
    }
    void update_settings(std::istream &is) {
        std::string identifier;
        // string scratch;
        while (!is.eof()) {
            is >> identifier;
            if (identifier.compare("hessian_shift_print_level:") == 0) {
                update_junk(is, map_hessian_shift_print_levels, hessian_shift_print_level);
            } else if (identifier.compare("initial_penalty_parameter:") == 0) {
                is >> initial_penalty_parameter;
            } else if (identifier.compare("qp_print_level:") == 0) {
                update_junk(is, map_qp_print_levels, qp_print_level);
            } else if (identifier.compare("termination_tolerance:") == 0) {
                is >> termination_tolerance;
            } else {
                std::cerr << "undefined identifier in input file: [" << identifier << "]." << std::endl;
            }
        }
    }
    
    void print() {
        std::cout << "All settings:" << std::endl
            << "hessian_shift_print_level_t hessian_shift_print_level = " << hessian_shift_print_level << std::endl
            << "qp_print_level_t qp_print_level = " << qp_print_level << std::endl
            << "double initial_penalty_parameter = " << initial_penalty_parameter << std::endl
            << "double termination_tolerance = " << termination_tolerance << std::endl
            << std::endl;
    }
    hessian_shift_print_level_t get_hessian_shift_print_level() { return hessian_shift_print_level; }
    double get_initial_penalty_parameter() { return initial_penalty_parameter; }
    qp_print_level_t get_qp_print_level() { return qp_print_level; }
    double get_termination_tolerance() { return termination_tolerance; }
    
private:
    hessian_shift_print_level_t hessian_shift_print_level;
    double initial_penalty_parameter;
    qp_print_level_t qp_print_level;
    double termination_tolerance;
    
    std::map<std::string, hessian_shift_print_level_t> map_hessian_shift_print_levels;
    std::map<std::string, qp_print_level_t> map_qp_print_levels;
};


#endif /* end of include guard: ISQO_CONFIG_CC_3N2PPVP3 */
