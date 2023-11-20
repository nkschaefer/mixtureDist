#ifndef _COMMON_H
#define _COMMON_H
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>

// ===== MISC. MATH FUNCTIONS =====

double factorial(int n);
double pois(int x, int l);
int binom_coef(int n, int k);
double dlgamma(double x, double a, double b);
double pgamma(double x, double a, double b);
std::pair<double, double> gamma_moments(double mu, double var);
double pois_cdf(int x, int l);
double expPDF(float x, float l);
double dbinom(double n, double k, double p);
double dlmultinom(const std::vector<double>& x, const std::vector<double>& p);
double dldirichlet(const std::vector<double>& x, const std::vector<double>& p);
double binom_cdf_lower(int n, int k, float p);
float pdfhyper(int N, int K, int n, int k);
float dnorm(float x, float mu, float sigma);
float cnorm(float x, float mu, float sigma);
float gaussPDF(float x, float mu, float sigma);
float betafunc(int a, int b);
float lbetaf(float a, float b);
long double betaf(float a, float b);
float dlbeta(float x, float a, float b);
float dbeta(float x, float a, float b);
float pbeta(float x, float a, float b);
std::pair<float, float> betafit(std::vector<float>& vals);
std::pair<float, float> betafit_disp(std::vector<float>& vals, float loc);
double binom_coef_log(double n, double k);
double dbetabin(double x, double n, double a, double b);
float dnbinom(int x, int mu, int phi);
float dnbinom_muvar(int x, float mu, float var);
float pnbinom(int x, int mu, int phi);
std::pair<int, float> nbinom_moments(float mean, float var);
std::pair<float, float> beta_loc_scale_to_a_b(float loc, float scale);
float binom_approx(int n, int k, float frac);
float binom_approx_scaled(int n, int k, float scale, float frac);
std::pair<double, double> welford(std::vector<float>& vals);
void gaussEM(std::vector<float>& obs, std::vector<float>& weights, std::vector<float>& mu,
    std::vector<float>& sigma, std::vector<short>& assignments);
float int2gauss(float mu1, float sigma1, float mu2, float sigma2);

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
        bool update( const std::vector<double>& means, const std::vector<double>& vars);

        // Visualize data in distribution
        void print();
        void print(int indentation_level);
        
        // How many parameters do all component distributions of this combined have?
        int get_n_parameters();

        // Create a new distribution, mapped to its name
        static void register_dist(std::string name, int n_inputs, int n_params, loglik_func ll, update_func u);
        
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
            int start_idx,
            int n_inputs,
            std::vector<double>& params,
            const std::vector<bool>& params_frozen,
            const bool all_frozen );
};

// MIXTURE MODELS =====

// Generalized function to handle globally shared parameters of component distributions.
// Runs just after all distributions have updated themselves (by fitting to the data).
// Gives pointers to all distributions in a set with shared parameters, and those 
// distributions' weights (scaled so that they sum to one within this set), so that an
// external program can reconcile these parameters based on some other, global 
// parameters they might depend on. It is up to that program to decide what to do and
// how, or if, to adjust each distribution's parameters.
typedef std::function< void( const std::vector<mixtureDist*>&,
    const std::vector<double>&,
    std::vector<double>& ) > shared_params_callback;

class mixtureModel{
    public:
        int n_components;
        std::vector<double> weights;
        std::vector<mixtureDist> dists;
        std::vector<short> assignments;
        float loglik;
        float bic;
        float aic;
        double** responsibility_matrix;
        
        // Constructor 
        mixtureModel();
        mixtureModel(std::vector<mixtureDist>& dists);
        mixtureModel(std::vector<mixtureDist>& dists, std::vector<double> weights);
        
        // Copy constructor
        mixtureModel(const mixtureModel& m);
        
        // Are there any global parameters that these shared distributions depend on? 
        // Store initial guesses in here and return to/update them in the callback function
        std::vector<std::vector<double> > shared_params;

        // Shortcuts for when using only one distribution type
        // Single-parameter distribution
        mixtureModel(std::string, std::vector<double> params);
        mixtureModel(std::string, std::vector<double> params, std::vector<double> weights);
        // Two-parameter distribution
        mixtureModel(std::string, std::vector<std::pair<double, double> > params);
        mixtureModel(std::string, std::vector<std::pair<double, double> > params, std::vector<double> weights);
        
        // Destructor
        ~mixtureModel();
        
        // Update mixture model parameters
        void set_delta_thresh(float delta_thresh);
        void set_maxits(int maxits);
        void set_verbosity(short level);
        
        // Set certain component distributions to share parameters
        bool set_shared_params(std::vector<int> shared_param_dists, shared_params_callback fun);
        bool set_shared_params(std::vector<int> shared_param_dists, shared_params_callback fun,
            std::vector<double> shared_meta_params);

        // Overloaded function to fit data passed in different ways
        double fit(const std::vector<std::vector<double> >& obs);
        double fit(const std::vector<std::vector<double> >& obs, std::vector<double>& obs_weights);
        
        double fit(const std::vector<std::vector<int> >& obs);
        double fit(const std::vector<std::vector<int> >& obs, std::vector<double>& obs_weights);

        double fit(const std::vector<double>& obs);
        double fit(const std::vector<double>& obs, std::vector<double>& obs_weights);

        double fit(const std::vector<int>& obs);
        double fit(const std::vector<int>& obs, std::vector<double>& obs_weights);

        void print();
    
    private:
        bool print_lls;
        bool print_dists;

        float delta_thresh;
        int maxits;
        int n_obs; // store size of current data set
        
        // Allow component distributions to share parameters. This will be
        // done by allowing each of these to update their parameters, 
        // and then computing a weighted average of each parameter, using
        // the distribution weights.
        std::vector<std::set<int> > shared_params_groups; 
        // What function should be called (externally) to reconcile distributions that share
        // parameters, after each maximization step?
        std::vector<shared_params_callback> shared_params_callbacks;
        // Set up a new mixture model
        void init(std::vector<mixtureDist>& dists, std::vector<double>& weights);

        // Allocate responsibility matrix for a data set
        void init_responsibility_matrix(int nobs);

        // De-allocate responsibility matrix
        void free_responsibility_matrix();
        
        // Compute log likelihood of a data set under current model
        double compute_loglik(const std::vector<std::vector<double> >& obs);

        // Run fitting for the case where there is only one distribution in the mixture model.
        double fit_single(const std::vector<std::vector<double> >& obs, std::vector<double>& obs_weights);

        // After fitting to data, compute the BIC and AIC
        void compute_bic(long int n_obs);
}; 

#endif
