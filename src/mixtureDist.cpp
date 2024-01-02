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
#include <functional>
#include <cstdlib>
#include <utility>
#include <math.h>
#include "functions.h"
#include "mixtureDist.h"

using std::cout;
using std::endl;
using namespace std;

// ===== mixtureDist =====

/**
 * Define log likelihood and update functions for pre-set probability distributions.
 */
double mixtureDist::ll_poisson(const vector<double>& input, 
    int start_idx,
    int n_inputs, 
    const vector<double>& params){
    return dpois(input[start_idx], params[0]);
}
bool mixtureDist::update_poisson(const vector<double>& means,
    const vector<double>& vars,
    int start_idx,
    int n_inputs,
    vector<double>& params,
    const vector<bool>& params_frozen,
    const bool all_frozen){
    if (means[start_idx] <= 0 || vars[start_idx] <= 0){
        return false;
    }
    if (!all_frozen && !params_frozen[0]){
        params[0] = means[start_idx];
    }
    return true;
}
double mixtureDist::ll_gamma(const vector<double>& input, 
    int start_idx,
    int n_inputs, 
    const vector<double>& params){
    return dgamma(input[start_idx], params[0], params[1]); 
}
bool mixtureDist::update_gamma(const vector<double>& means,
    const vector<double>& vars,
    int start_idx,
    int n_inputs,
    vector<double>& params,
    const vector<bool>& params_frozen,
    const bool all_frozen){
    if (means[start_idx] <= 0 || vars[start_idx] <= 0){
        return false;
    }    
    pair<double, double> newparams = gamma_moments(means[start_idx], vars[start_idx]);
    if (isinf(newparams.first) || isinf(newparams.second) || newparams.first > 1e8 || newparams.second > 1e8){
        return false;
    }
    if (newparams.first <= 0 || newparams.second <= 0){
        return false;
    }
    if (!all_frozen && !params_frozen[0]){
        params[0] = newparams.first;
    }
    if (!all_frozen && !params_frozen[1]){
        params[1] = newparams.second;
    }
    return true;
}
double mixtureDist::ll_beta(const vector<double>& input, 
    int start_idx,
    int n_inputs, 
    const vector<double>& params){
    return dlbeta(input[start_idx], params[0], params[1]);
}
bool mixtureDist::update_beta(const vector<double>& means,
    const vector<double>& vars,
    int start_idx,
    int n_inputs,
    vector<double>& params,
    const vector<bool>& params_frozen,
    const bool all_frozen){
    if (means[start_idx] <= 0 || means[start_idx] == 1){
        return false;
    }
    if (vars[start_idx] == 0){
        return false;
    }
    double mean = means[start_idx];
    double var = vars[start_idx];
    double alpha = mean*((mean*(1-mean))/var - 1);
    double beta = (alpha * (1-mean))/mean; 
    if (isinf(alpha) || alpha > 1e6 || isinf(beta) || beta > 1e6){
        return false;
    }
    if (!params_frozen[0] && !all_frozen){
        params[0] = alpha;
    }
    if (!params_frozen[1] && !all_frozen){
        params[1] = beta;
    }
    return true; 
}
double mixtureDist::ll_gauss(const vector<double>& input, 
    int start_idx,
    int n_inputs, 
    const vector<double>& params){
    return dnorm(input[start_idx], params[0], params[1]); 
}

bool mixtureDist::update_gauss(const vector<double>& means,
    const vector<double>& vars,
    int start_idx,
    int n_inputs,
    vector<double>& params,
    const vector<bool>& params_frozen,
    const bool all_frozen){

    if (vars[start_idx] == 0){
        // Normal distribution with zero variance is illegal
        return false;
    }
    if (!all_frozen && !params_frozen[0]){
        params[0] = means[start_idx];
    }
    if (!all_frozen && !params_frozen[1]){
        params[1] = sqrt(vars[start_idx]);
    }
    return true;
}

double mixtureDist::ll_binom(const vector<double>& input, 
    int start_idx,
    int n_inputs, 
    const vector<double>& params){
    // input 1 is n (successes), input 2 is k (total)
    return dbinom(input[start_idx+1], input[start_idx], params[0]); 
}

bool mixtureDist::update_binom(const vector<double>& means,
    const vector<double>& vars,
    int start_idx,
    int n_inputs,
    vector<double>& params,
    const vector<bool>& params_frozen,
    const bool all_frozen){
    
    // This distribution takes 2 inputs, so we will consume two of the means
    // to compute the expectation.
    double mean1 = means[start_idx];
    double mean2 = means[start_idx+1];
    if (mean1 == 0 || mean1 == mean2){
        // This will result in expectation of 0 or 1, which is impossible
        return false;
    }
    else if (!all_frozen && !params_frozen[0]){
        params[0] = mean1/mean2;
    }
    return true;
}

double mixtureDist::ll_multinom(const vector<double>& input, 
    int start_idx,
    int n_inputs, 
    const vector<double>& params){
    
    // Build a vector of all input values
    vector<double> inputs_this;
    // Ensure parameters sum to 1
    double param_tot = 0.0;
    for (int i = 0; i < n_inputs; ++i){
        inputs_this.push_back(input[start_idx + i]);
        param_tot += params[i];
    }
    //if (param_tot != 1.0){
    //    return 1.0/0.0;
    //}
    // Params should already contain the correct number of parameters.
    if (params.size() != n_inputs){
        // Return an invalid value so we know something went wrong.
        return 1.0/0.0;
    }
    return dmultinom(inputs_this, params);
}

bool mixtureDist::update_multinom(const vector<double>& means,
    const vector<double>& vars,
    int start_idx,
    int n_inputs,
    vector<double>& params,
    const vector<bool>& params_frozen,
    const bool all_frozen){
    
    if (all_frozen){
        return true;
    }

    double tot_p = 0.0;
    double tot = 0.0;
    for (int i = 0; i < n_inputs; ++i){
        // Leave frozen parameters out of the sum
        if (!params_frozen[i]){
            tot += means[start_idx+i];
            tot_p += params[i];
        }
    }
    bool all_valid = true;
    for (int i = 0; i < n_inputs; ++i){
        if (!params_frozen[i]){
            params[i] = (means[start_idx+i]/tot)*tot_p;
        }    
        if (params[i] == 0.0 || params[i] == 1.0){
            all_valid = false;
        }
    }
    return all_valid;
}

std::map<std::string, int> mixtureDist::registered_n_inputs;
std::map<std::string, int> mixtureDist::registered_n_params;
std::map<std::string, loglik_func> mixtureDist::registered_ll_func;
std::map<std::string, update_func> mixtureDist::registered_update_func;

// Register all preset distributions with class
void mixtureDist::register_dist(string name, int n_inputs, int n_params, loglik_func llf, update_func uf){
    mixtureDist::registered_n_inputs.insert(make_pair(name, n_inputs));
    mixtureDist::registered_n_params.insert(make_pair(name, n_params));
    mixtureDist::registered_ll_func.insert(make_pair(name, llf));
    mixtureDist::registered_update_func.insert(make_pair(name, uf));
}

/**
 * Set up pre-set distributions, whenever a new class member is created
 */
void mixtureDist::auto_register(){
    if (mixtureDist::registered_n_inputs.count("poisson") == 0){    
        mixtureDist::register_dist("poisson", 1, 1, mixtureDist::ll_poisson, mixtureDist::update_poisson);
    }
    if (mixtureDist::registered_n_inputs.count("gamma") == 0){
        mixtureDist::register_dist("gamma", 1, 2, mixtureDist::ll_gamma, mixtureDist::update_gamma);
    }
    if (mixtureDist::registered_n_inputs.count("beta") == 0){
        mixtureDist::register_dist("beta", 1, 2, mixtureDist::ll_beta, mixtureDist::update_beta);
    }
    if (mixtureDist::registered_n_inputs.count("gauss") == 0){
        mixtureDist::register_dist("gauss", 1, 2, mixtureDist::ll_gauss, mixtureDist::update_gauss);
    }
    if (mixtureDist::registered_n_inputs.count("normal") == 0){
        mixtureDist::register_dist("normal", 1, 2, mixtureDist::ll_gauss, mixtureDist::update_gauss);
    }
    if (mixtureDist::registered_n_inputs.count("binomial") == 0){
        mixtureDist::register_dist("binomial", 2, 1, mixtureDist::ll_binom, mixtureDist::update_binom);
    }
    if (mixtureDist::registered_n_inputs.count("multinomial") == 0){
        mixtureDist::register_dist("multinomial", -1, -1, mixtureDist::ll_multinom, mixtureDist::update_multinom);
    }
}

mixtureDist::mixtureDist(){
    this->n_components = 0;
    this->frozen = false;
    this->name = "";
    mixtureDist::auto_register();
}

mixtureDist::~mixtureDist(){
    
}

mixtureDist::mixtureDist(const mixtureDist& m){
    // Should not be needed, since the other object was already created
    //mixtureDist::auto_register();

    this->names.clear();
    this->weights.clear();
    this->n_inputs.clear();
    this->loglik_funcs.clear();
    this->update_funcs.clear();
    this->n_params.clear();
    this->params.clear();
    this->params_frozen.clear();
    
    this->n_components = m.n_components;
    this->frozen = m.frozen;
    this->name = m.name;

    for (int i = 0; i < m.n_components; ++i){
        this->names.push_back(m.names[i]);
        this->weights.push_back(m.weights[i]);
        this->n_inputs.push_back(m.n_inputs[i]);
        this->loglik_funcs.push_back(m.loglik_funcs[i]);
        this->update_funcs.push_back(m.update_funcs[i]);
        this->n_params.push_back(m.n_params[i]);
        this->params.push_back(m.params[i]);
        this->params_frozen.push_back(m.params_frozen[i]);
    }
}

void mixtureDist::print_unregistered(string name){
    fprintf(stderr, "ERROR: %s is not a recognized distribution.\n", name.c_str());
    fprintf(stderr, "Please choose a pre-registered distribution or create one yourself.\n");
    fprintf(stderr, "Registered distributions:\n");
    for (map<string, int>::iterator r = mixtureDist::registered_n_inputs.begin(); r != 
        mixtureDist::registered_n_inputs.end(); ++r){
        fprintf(stderr, "\t%s\n", r->first.c_str());
    }
    fprintf(stderr, "To register your own, call mixtureDist::register_dist() with the following arguments:\n");
    fprintf(stderr, "\tname (string)\n");
    fprintf(stderr, "\tnumber of inputs (int), e.g. 1 for univariate distributions\n");
    fprintf(stderr, "\tnumber of parameters (int)\n");
    fprintf(stderr, "\tlikelihood function of the form:\n");
    fprintf(stderr, "\t\tdouble func(const vector<double>& input, int index_in_input, const vector<double>& parameters)\n");
    fprintf(stderr, "\tfunction to update parameters from data moments of the form:\n");
    fprintf(stderr, "\t\tbool func(const vector<double>& means, const vector<double>& variances, int index_in_means_and_variances, vector<double>& parameters\n");
    fprintf(stderr, "\t\tthis function should return true on success and false on failure.\n");
}

bool mixtureDist::add_dist(string name,
    vector<double> params_input, 
    double weight, 
    vector<bool> params_input_frozen){

    if (mixtureDist::registered_n_inputs.count(name) == 0){
        mixtureDist::print_unregistered(name);
        return false;
    }
    else{

        this->n_components++;
        this->n_inputs.push_back(mixtureDist::registered_n_inputs[name]);
        this->names.push_back(name);
        this->n_params.push_back(mixtureDist::registered_n_params[name]);
        // Check that correct number of parameters was input
        if (mixtureDist::registered_n_params[name] == -1){
            // warn user to set the number of inputs (flexible distribution)?
            //fprintf(stderr, "WARNING: remember to set n_inputs\n");
        }
        else if (params_input.size() != mixtureDist::registered_n_params[name]){
            fprintf(stderr, "ERROR: %ld parameters given for %s distribution, which requires \
%d parameters.\n", params_input.size(), name.c_str(), mixtureDist::registered_n_params[name]);
            return false;
        }
        this->params.push_back(params_input);
        this->weights.push_back(weight);
        if (params_input_frozen.size() != params_input.size()){
            fprintf(stderr, "ERROR: %ld parameters passed, but frozen settings for %ld parameters given.\n",
                params_input.size(), params_input_frozen.size());
            return false;
        }
        this->params_frozen.push_back(params_input_frozen);
        this->loglik_funcs.push_back(mixtureDist::registered_ll_func[name]);
        this->update_funcs.push_back(mixtureDist::registered_update_func[name]);
        return true;
    }
}

bool mixtureDist::add_dist(string name, vector<double> params_input, double weight){
    vector<bool> params_input_frozen;
    for (int i = 0; i < params_input.size(); ++i){
        params_input_frozen.push_back(false);
    }
    return this->add_dist(name, params_input, weight, params_input_frozen);
}

bool mixtureDist::add_dist(string name, vector<double> params_input, vector<bool> params_input_frozen){
    return this->add_dist(name, params_input, 1.0, params_input_frozen);
}

bool mixtureDist::add_dist(string name, vector<double> params_input){
    vector<bool> params_frozen;
    for (int i = 0; i < params_input.size(); ++i){
        params_frozen.push_back(false);
    }
    return this->add_dist(name, params_input, 1.0, params_frozen);
}

bool mixtureDist::add_dist(string name, double weight){
    vector<double> params;
    vector<bool> params_frozen;
    if (mixtureDist::registered_n_params.count(name) > 0){
        for (int i = 0; i < mixtureDist::registered_n_params[name]; ++i){
            params.push_back(0.0);
            params_frozen.push_back(false);
        }
    }
    return this->add_dist(name, params, weight, params_frozen);
}

bool mixtureDist::add_dist(string name){
    // Initialize parameters to zeroes as placeholders
    vector<double> params;
    if (mixtureDist::registered_n_params.count(name) > 0){
        for (int i = 0; i < mixtureDist::registered_n_params[name]; ++i){
            params.push_back(0.0);
        }
    }
    return this->add_dist(name, params, 1.0);
}

mixtureDist::mixtureDist(string name){
    this->n_components = 0;
    mixtureDist::auto_register();
    if (! this->add_dist(name)){
        exit(1);
    }
    this->frozen = false;
}

mixtureDist::mixtureDist(string name, vector<double> params_input){
    this->n_components = 0;
    mixtureDist::auto_register();
    if (!this->add_dist(name, params_input)){
        exit(1);
    }
    this->frozen = false;
}

/**
 * Shortcut to initialize single-parameter distributions
 */
mixtureDist::mixtureDist(string name, double param){
    this->n_components = 0;
    mixtureDist::auto_register();
    vector<double> params_input;
    params_input.push_back(param);
    if (!this->add_dist(name, params_input)){
        exit(1);
    }
    this->frozen = false;
}

/**
 * Shortcut to initialize single-parameter distributions with frozen setting
 */
mixtureDist::mixtureDist(string name, double param, bool frozen){
    this->n_components = 0;
    mixtureDist::auto_register();
    vector<double> params_input;
    params_input.push_back(param);
    vector<bool> params_input_frozen;
    params_input_frozen.push_back(frozen);
    if (!this->add_dist(name, params_input, params_input_frozen)){
        exit(1);
    }
    this->frozen = false;
}

/**
 * Shortcut to initialize two-parameter distributions
 */
mixtureDist::mixtureDist(string name, double param1, double param2){
    this->n_components = 0;
    mixtureDist::auto_register();
    vector<double> params_input;
    params_input.push_back(param1);
    params_input.push_back(param2);
    if (!this->add_dist(name, params_input)){
        exit(1);
    }
    this->frozen = false;
}

/**
 * Shortcut to initialize two-parameter distributions with frozen settings
 */
mixtureDist::mixtureDist(string name, double param1, double param2, bool frozen1, bool frozen2){
    this->n_components = 0;
    mixtureDist::auto_register();
    vector<double> params_input;
    params_input.push_back(param1);
    params_input.push_back(param2);
    vector<bool> params_input_frozen;
    params_input_frozen.push_back(frozen1);
    params_input_frozen.push_back(frozen2);
    if (!this->add_dist(name, params_input, params_input_frozen)){
        exit(1);
    }
    this->frozen = false;
}
// Initialize a compound distribution

mixtureDist::mixtureDist(vector<string> names){
    this->n_components = 0;
    mixtureDist::auto_register();
    for (int i = 0; i < names.size(); ++i){
        if (!this->add_dist(names[i])){
            exit(1);
        }
    }
    this->frozen = false;
}

mixtureDist::mixtureDist(vector<string> names, vector<vector<double> > params_input){
    this->n_components = 0;
    mixtureDist::auto_register();
    for (int i = 0; i < names.size(); ++i){
        if (!this->add_dist(names[i], params_input[i])){
            exit(1);
        }
    }
    this->frozen = false;
}

mixtureDist::mixtureDist(vector<string> names, vector<vector<double> > params_input,
    vector<vector<bool> > params_input_frozen){
    this->n_components = 0;
    mixtureDist::auto_register();
    for (int i = 0; i < names.size(); ++i){
        if (!this->add_dist(names[i], params_input[i], params_input_frozen[i])){
            exit(1);
        }
    }
    this->frozen = false;
}

mixtureDist::mixtureDist(vector<string> names, vector<double> weights){
    this->n_components = 0;
    mixtureDist::auto_register();
    for (int i = 0; i < names.size(); ++i){
        if (!this->add_dist(names[i], weights[i])){
            exit(1);
        }
    }
    this->frozen = false;
}

mixtureDist::mixtureDist(vector<string> names, vector<vector<double> > params_input, vector<double> weights){
    this->n_components = 0;
    mixtureDist::auto_register();
    for (int i = 0; i < names.size(); ++i){
        if (!this->add_dist(names[i], params_input[i], weights[i])){
            exit(1);
        }
    }
    this->frozen = false;
}

mixtureDist::mixtureDist(vector<string> names, vector<vector<double> > params_input,
    vector<double> weights, vector<vector<bool> > params_input_frozen){
    this->n_components = 0;
    mixtureDist::auto_register();
    for (int i = 0; i < names.size(); ++i){
        if (!this->add_dist(names[i], params_input[i], weights[i], params_input_frozen[i])){
            exit(1);
        }
    }
    this->frozen = false;
}

double mixtureDist::loglik( const vector<double>& inputs ){
    double ll = 0.0;
    int start_idx = 0;
    for (int i = 0; i < this->n_components; ++i){
        if (this->n_inputs[i] == -1 ){
           fprintf(stderr, "ERROR: number of inputs not set for a distribution. \
If using a distribution with a flexible number of inputs, you must call set_num_inputs() \
to tell it how many values to take in.\n");
           this->print();
           exit(1);
           
        }
        ll += log2(this->weights[i]) + 
            this->loglik_funcs[i]( inputs, start_idx, this->n_inputs[i], this->params[i] );
            start_idx += this->n_inputs[i];
    }
    return ll;
}

bool mixtureDist::update( const vector<double>& means, const vector<double>& vars ){
    bool success = true;
    int start_idx = 0;
    for (int i = 0; i < this->n_components; ++i){
        if (this->n_inputs[i] == -1 ){
            fprintf(stderr, "ERROR: number of inputs not set for a distribution. \
If using a distribution with a flexible number of inputs, you must call set_num_inputs() \
to tell it how many values to take in.\n");
            exit(1);
        }
        success = success & this->update_funcs[i]( means, vars, start_idx, this->n_inputs[i],
            this->params[i], this->params_frozen[i], this->frozen );
        start_idx += this->n_inputs[i];  
    }
    return success;
} 

void mixtureDist::print_one(string& indent_str, int idx, bool incl_weight){
    string indent_str2 = indent_str + "   ";
    if (this->n_components <= 1 && this->name != ""){
        indent_str += this->name + ": ";
    }
    if (incl_weight){
        fprintf(stderr, "%s%s distribution weight = %.3f\n", indent_str.c_str(),
            this->names[idx].c_str(), this->weights[idx]);
    }
    else{
        fprintf(stderr, "%s%s distribution\n", indent_str.c_str(),
            this->names[idx].c_str());
    }
    fprintf(stderr, "%sParameters: [", indent_str2.c_str());
    int npar = this->n_params[idx];
    if (npar == -1){
        npar = this->n_inputs[idx];
    }
    for (int i = 0; i < npar; ++i){
        fprintf(stderr, " %.3f", this->params[idx][i]);
        if (this->params_frozen[idx][i]){
            fprintf(stderr, "*");
        }
    }
    fprintf(stderr, " ]\n");
}

void mixtureDist::print(int indentation_level){
    string indent_str = "";
    for (int i = 0; i < indentation_level; ++i){
        indent_str += "   ";
    }
    if (this->n_components > 1){
        string namestr = "";
        if (this->name != ""){
            namestr = this->name + ": ";
        }
        string frzstr = "";
        if (this->frozen){
            frzstr = "frozen ";
        }
        fprintf(stderr, "%s%s%sDist with %d components:\n", frzstr.c_str(), 
            indent_str.c_str(), namestr.c_str(), this->n_components);
        indent_str += "   ";
        // Don't print weights if they're identical
        double weight1 = this->weights[0];
        bool weights_match = true;
        for (int i = 1; i < this->n_components; ++i){
            if (this->weights[i] != weight1){
                weights_match = false;
                break;
            }
        }
        for (int i = 0; i < this->n_components; ++i){
            this->print_one(indent_str, i, !weights_match);
        }
    }
    else{
        this->print_one(indent_str, 0, false);
    }
}

void mixtureDist::print(){
    this->print(0); 
}

int mixtureDist::get_n_parameters(){
    int nparam = 0;
    for (int i = 0; i < this->n_components; ++i){
        nparam += this->n_params[i];
    }
    // component weights are not a free/fitted parameter, so ignore those
    return nparam;
}

void mixtureDist::set_num_inputs(int component_idx, int ni){
    if (component_idx > this->n_components-1){
        fprintf(stderr, "ERROR: cannot set number of inputs on component %d; \
distribution only has %d components\n", component_idx, this->n_components);
        exit(1);
    }
    else if (this->n_inputs[component_idx] != -1 && this->n_inputs[component_idx] != 
        ni){
        fprintf(stderr, "ERROR: component %d is set to take %d inputs already. \
Can only change number of inputs for flexible distributions with no default number.\n",
            component_idx, ni);
        exit(1);
    }
    this->n_inputs[component_idx] = ni;
}

void mixtureDist::set_num_inputs(int ni){
    // Assume this is single-component; use as a shortcut for setting num inputs for 
    // single component.
    this->set_num_inputs(0, ni);
}

