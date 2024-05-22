#ifndef _MIXTUREDIST_H
#define _MIXTUREDIST_H
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <functional>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>

// ===== COMPONENT DISTRIBUTIONS FOR MIXTURE MODELS =====

// Generalized function to compute log likelihood of input data (1+ dimensional) under
// a probability distribution
typedef std::function< double( const std::vector<double>&, 
    int,
    int, 
    const std::vector<double>& ) > loglik_func;

// Generalized function to update the parameters of a probability distribution using
// method of moments, given a vector of means and variances (1+ dimensions)
typedef std::function< bool( const std::vector<double>&, 
    const std::vector<double>&, 
    const std::map<std::pair<int, int>, double>&,
    int,
    int, 
    std::vector<double>&,
    const std::vector<bool>&,
    const bool ) > update_func;

class mixtureDist{
    public: 
        
        // Optional: name (will be printed and can be accessed by the user)
        std::string name;

        // How many parameters for each component distribution?
        std::vector<int> n_params;
        
        // All parameters for each component distribution
        std::vector<std::vector<double> > params;
       
        // Lock all parameters associated with this distribution
        bool frozen;
        
        // Needs to know about covariance between inputs?
        bool needs_cov;

        // Allow any or all parameters of this distribution to be frozen 
        // (will not update to fit data)
        std::vector<std::vector<bool> > params_frozen;

        // Initialize with no parameters 
        mixtureDist();
        
        // Destructor
        ~mixtureDist();
        
        // Copy constructor
        mixtureDist(const mixtureDist& m);

        // Initialize for a single distribution
        // Weights are nonsensical when not a compound distribution - 
        // don't pass weights and set all to 1
        mixtureDist(std::string name);
        mixtureDist(std::string name, double param);
        mixtureDist(std::string name, double param, bool param_frozen);
        mixtureDist(std::string name, double param1, double param2);
        mixtureDist(std::string name, double param1, double param2, bool param1_frozen, bool param2_frozen);
        mixtureDist(std::string name, std::vector<double> params);
        mixtureDist(std::string name, std::vector<double> params, std::vector<bool> params_frozen);

        // Initialize for compound, independent distributions
        mixtureDist(std::vector<std::string> names);
        mixtureDist(std::vector<std::string> names, std::vector<std::vector<double> > params);
        mixtureDist(std::vector<std::string> names, std::vector<std::vector<double> > params,
            std::vector<std::vector<bool> > params_input_frozen);
        mixtureDist(std::vector<std::string> names, std::vector<double> weights);
        mixtureDist(std::vector<std::string> names, std::vector<std::vector<double> > params, 
            std::vector<double> weights);
        mixtureDist(std::vector<std::string> names, std::vector<std::vector<double> > params,
            std::vector<double> weights, std::vector<std::vector<bool> > params_input_frozen);

        // Add another distribution
        bool add_dist(std::string name);
        bool add_dist(std::string name, std::vector<double> params);
        bool add_dist(std::string name, double weight);
        bool add_dist(std::string name, std::vector<double> params, double weight);
        bool add_dist(std::string name, std::vector<double> params, std::vector<bool> params_frozen);
        bool add_dist(std::string name, std::vector<double> params, double weight, 
            std::vector<bool> params_frozen); 
        
        // Get log likelihood of one observation (or multiple, if multivariate or compound distribution)
        double loglik( const std::vector<double>& );
        // Given summaries of the data, update parameters of the distribution using MoM
        bool update( const std::vector<double>& means, const std::vector<double>& vars,
            const std::map<std::pair<int, int>, double>& covs);

        // Visualize data in distribution
        void print();
        void print(int indentation_level);
        
        // How many parameters do all component distributions of this combined have?
        int get_n_parameters();

        // Create a new distribution, mapped to its name
        static void register_dist(std::string name, int n_inputs, int n_params, bool nc,
            loglik_func ll, update_func u);
        
        // If using a flexible distribution (such as multinomial), which can handle a non-fixed
        // number of inputs, tell the distribution here how many inputs to expect. 
        // Sets n_i inputs for component index comp_idx.
        void set_num_inputs(int comp_idx, int n_i);    
        // Above function but assume one component
        void set_num_inputs(int n_i);

    private:
        // How many sub-distributions (> 1 if compound distribution)
        int n_components;
        
        // Store names of component distributions
        std::vector<std::string> names;

        // How many inputs for each component distribution?
        std::vector<int> n_inputs;
        
        // How is each component distribution weighted?
        std::vector<double> weights;
        
        // What is the actual function used to calculate likelihood under each component?
        std::vector<loglik_func> loglik_funcs;

        // What is the actual function used to update parameters of each component to fit data?
        std::vector<update_func> update_funcs;
        
        // Keep a map of default function names to information about them
        // These will be accessed using the register_dist() static class method.
        static std::map<std::string, int> registered_n_inputs;
        static std::map<std::string, int> registered_n_params;
        static std::map<std::string, bool> registered_needs_cov;
        static std::map<std::string, loglik_func> registered_ll_func;
        static std::map<std::string, update_func> registered_update_func;
        static void auto_register();

        static void print_unregistered(std::string name);
        void print_one(std::string& indent_str, int idx, bool print_weight);

        static double ll_poisson( const std::vector<double>& input, 
            int start_idx, 
            int n_inputs,
            const std::vector<double>& params );
        static bool update_poisson( const std::vector<double>& means, 
            const std::vector<double>& vars,
            const std::map<std::pair<int, int>, double>& covs,
            int start_idx,
            int n_inputs,
            std::vector<double>& params,
            const std::vector<bool>& params_frozen,
            const bool all_frozen );
        static double ll_gamma( const std::vector<double>& input, 
            int start_idx,
            int n_inputs, 
            const std::vector<double>& params );
        static bool update_gamma( const std::vector<double>& means,
            const std::vector<double>& vars,
            const std::map<std::pair<int, int>, double>& covs,
            int start_idx,
            int n_inputs,
            std::vector<double>& params,
            const std::vector<bool>& params_frozen,
            const bool all_frozen );
        static double ll_beta( const std::vector<double>& input, 
            int start_idx,
            int n_inputs, 
            const std::vector<double>& params );
        static bool update_beta( const std::vector<double>& means, 
            const std::vector<double>& vars,
            const std::map<std::pair<int, int>, double>& covs,
            int start_idx,
            int n_inputs,
            std::vector<double>& params,
            const std::vector<bool>& params_frozen,
            const bool all_frozen );
        static double ll_gauss( const std::vector<double>& input, 
            int start_idx,
            int n_inputs, 
            const std::vector<double>& params );
        static bool update_gauss( const std::vector<double>& means, 
            const std::vector<double>& vars,
            const std::map<std::pair<int, int>, double>& covs,
            int start_idx,
            int n_inputs,
            std::vector<double>& params,
            const std::vector<bool>& params_frozen,
            const bool all_frozen );
        static double ll_binom( const std::vector<double>& input,
            int start_idx,
            int n_inputs,
            const std::vector<double>& vars);
        static bool update_binom( const std::vector<double>& means,
            const std::vector<double>& vars,
            const std::map<std::pair<int, int>, double>& covs,
            int start_idx,
            int n_inputs,
            std::vector<double>& params,
            const std::vector<bool>& params_frozen,
            const bool all_frozen );
        static double ll_multinom( const std::vector<double>& input,
            int start_idx,
            int n_inputs,
            const std::vector<double>& vars);
        static bool update_multinom( const std::vector<double>& means,
            const std::vector<double>& vars,
            const std::map<std::pair<int, int>, double>& covs,
            int start_idx,
            int n_inputs,
            std::vector<double>& params,
            const std::vector<bool>& params_frozen,
            const bool all_frozen );
        static double ll_2dgauss( const std::vector<double>& input,
            int start_idx,
            int n_inputs,
            const std::vector<double>& vars);
        static bool update_2dgauss( const std::vector<double>& means,
            const std::vector<double>& vars,
            const std::map<std::pair<int, int>, double>& covs,
            int start_idx,
            int n_inputs,
            std::vector<double>& params,
            const std::vector<bool>& params_frozen,
            const bool all_frozen );
};

#endif
