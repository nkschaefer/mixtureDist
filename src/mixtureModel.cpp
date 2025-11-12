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
#include "mixtureDist.h"
#include "mixtureModel.h"

using std::cout;
using std::endl;
using namespace std;

// ===== mixtureModel =====

void mixtureModel::init(vector<mixtureDist>& dists_init,
    vector<double>& weights){
    
    this->n_components = dists_init.size(); 
    this->dists = dists_init; 
    this->aic = 0.0;
    this->bic = 0.0;
    this->loglik = 0.0;
    this->weights = weights;
    
    this->delta_thresh = 0.1;
    this->maxits = 1000;
    
    this->n_obs = -1;
    this->responsibility_matrix = NULL;
    
    this->print_lls = false;
    this->print_dists = false;
    
    this->has_callback_fun = false;
    this->e_only = false;
}

void mixtureModel::init_responsibility_matrix(int n_obs){
    if (this->n_obs != -1){
        // Responsibility matrix was formerly initialized for a different
        // data set
        this->free_responsibility_matrix();
    }
    
    this->n_obs = n_obs;

    // Store the weight of each (component, observation) combination
    //this->responsibility_matrix = (double*)malloc(sizeof(double[n_obs][n_components]));
    
    this->responsibility_matrix = new double*[this->n_obs];
    for (int i = 0; i < this->n_obs; ++i){
        this->responsibility_matrix[i] = new double[this->n_components];
    }
    
}

void mixtureModel::free_responsibility_matrix(){
    if (this->n_obs != -1){
        for (int i = 0; i < this->n_obs; ++i){
            delete[] this->responsibility_matrix[i];
            this->responsibility_matrix[i] = NULL;
        }
        delete[] this->responsibility_matrix;
        this->responsibility_matrix = NULL;
    }
    this->n_obs = -1;
}

mixtureModel::mixtureModel(){
    this->n_components = -1;
    this->responsibility_matrix = NULL;
    this->n_obs = -1;
    this->aic = 0.0;
    this->bic = 0.0;
    this->loglik = 0.0;
    this->delta_thresh = 0.1;
    this->maxits = 1000;
    this->print_lls = false;
    this->print_dists = false;
    this->has_callback_fun = false;
    this->e_only = false;
}

mixtureModel::mixtureModel(const mixtureModel& m){
    this->n_components = m.n_components;
    //this->free_responsibility_matrix();
    this->responsibility_matrix = NULL;
    this->n_obs = -1;
    if (m.n_obs != -1 && m.responsibility_matrix != NULL){
        this->init_responsibility_matrix(m.n_obs);
        for (int i = 0; i < m.n_obs; ++i){
            for (int j = 0; j < m.n_components; ++j){
                this->responsibility_matrix[i][j] = m.responsibility_matrix[i][j];
            }
        }
    }
    this->loglik = m.loglik;
    this->bic = m.bic;
    this->aic = m.aic;
    this->delta_thresh = m.delta_thresh;
    this->maxits = m.maxits;
    this->assignments = m.assignments;
    
    this->callback_fun = m.callback_fun;
    this->has_callback_fun = m.has_callback_fun; 
    this->shared_params = m.shared_params;

    for (int i = 0; i < m.n_components; ++i){
        this->weights.push_back(m.weights[i]);
        this->dists.push_back(m.dists[i]);
    }
    this->print_lls = false;
    this->print_dists = false;
    this->e_only = m.e_only;
}

mixtureModel::mixtureModel(vector<mixtureDist>& dists, vector<double> weights){
    if (dists.size() != weights.size()){
        fprintf(stderr, "ERROR: distribution and weight vectors must be the same size.\n");
        exit(1);
    }
    // Ensure weights sum to one
    double wsum = 0.0;
    for (int i = 0; i < weights.size(); ++i){
        wsum += weights[i];
    }
    if (wsum != 1.0){
        for (int i = 0; i < weights.size(); ++i){
            weights[i] /= wsum;
        }
    }
    this->init(dists, weights);
}

mixtureModel::mixtureModel(vector<mixtureDist>& dists){
    // Assign equal weight to each component
    vector<double> weights;
    for (int i = 0; i < dists.size(); ++i){
        weights.push_back(1.0 / (double)dists.size());
    }
    this->init(dists, weights);
}

mixtureModel::mixtureModel(string dist_name, vector<double> param_single){
    vector<double> weights;
    vector<mixtureDist> dists;
    for (int i = 0; i < param_single.size(); ++i){
        weights.push_back(1.0 / (double)param_single.size());
        dists.push_back(mixtureDist(dist_name, param_single[i]));
    }
    this->init(dists, weights); 
}

mixtureModel::mixtureModel(string dist_name, vector<double> param_single, vector<double> weights){
    if (param_single.size() != weights.size()){
        fprintf(stderr, "ERROR: parameter and weight vectors must be the same size.\n");
        exit(1);
    }
    double wsum = 0.0;
    vector<mixtureDist> dists;
    for (int i = 0; i < param_single.size(); ++i){
        wsum += weights[i];
        dists.push_back(mixtureDist(dist_name, param_single[i]));
    }
    for (int i = 0; i < param_single.size(); ++i){
        weights[i] /= wsum;
    }
    this->init(dists, weights);
}

mixtureModel::mixtureModel(string dist_name, vector<pair<double, double> > param_double){
    vector<double> weights;
    vector<mixtureDist> dists;
    for (int i = 0; i < param_double.size(); ++i){
        weights.push_back(1.0/(double)param_double.size());
        dists.push_back(mixtureDist(dist_name, param_double[i].first, param_double[i].second));
    }
    this->init(dists, weights);
}

mixtureModel::mixtureModel(string dist_name, vector<pair<double, double> > param_double, vector<double> weights){
    if (param_double.size() != weights.size()){
        fprintf(stderr, "ERROR: parameter vector and weight vector must be the same size.\n");
        exit(1);
    }
    double wsum = 0.0;
    vector<mixtureDist> dists;
    for (int i = 0; i < param_double.size(); ++i){
        wsum += weights[i];
        dists.push_back(mixtureDist(dist_name, param_double[i].first, param_double[i].second));
    }
    for (int i = 0; i < param_double.size(); ++i){
        weights[i] /= wsum;
    }
    this->init(dists, weights);
}

/**
 * Destructor
 */
mixtureModel::~mixtureModel(){
    this->weights.clear();
    this->assignments.clear();
    this->free_responsibility_matrix();
}

void mixtureModel::set_delta_thresh(float delta_thresh){
    this->delta_thresh = delta_thresh;
}

void mixtureModel::set_maxits(int maxits){
    this->maxits = maxits;
}

void mixtureModel::set_verbosity(short level){
    if (level == 0){
        this->print_lls = false;
        this->print_dists = false;
    }
    else if (level == 1){
        this->print_lls = true;
        this->print_dists = false;
    }
    else if (level == 2){
        this->print_lls = true;
        this->print_dists = true;
    }
    else{
        fprintf(stderr, "ERROR: invalid verbosity level %d\n", level);
        exit(1);
    }
}

void mixtureModel::freeze_dists(){
    this->e_only = true;
}

bool mixtureModel::set_callback(callback callback_func){
    vector<double> meta_params_empty;
    return this->set_callback(callback_func, meta_params_empty);
}

bool mixtureModel::set_callback(callback callback_func,
    vector<double> meta_params){
    
    this->has_callback_fun = true;
    this->callback_fun = callback_func;
    this->shared_params = meta_params;
    return true;
}

double mixtureModel::fit(const vector<double>& obs){
    vector<vector<double> > obs_multi;
    vector<double> obs_weights;
    for (int i = 0; i < obs.size(); ++i){
        vector<double> v;
        v.push_back(obs[i]);
        obs_multi.push_back(v);
        obs_weights.push_back(1.0);
    }
    return this->fit(obs_multi, obs_weights);
}

double mixtureModel::fit(const vector<double>& obs, vector<double>& obs_weights){
    vector<vector<double> > obs_multi;
    for (int i = 0; i < obs.size(); ++i){
        vector<double> v;
        v.push_back(obs[i]);
        obs_multi.push_back(v);
    }
    return this->fit(obs_multi, obs_weights);
}

double mixtureModel::fit(const vector<int>& obs){
    vector<vector<double> > obs_multi;
    vector<double> obs_weights;
    for (int i = 0; i < obs.size(); ++i){
        vector<double> v;
        v.push_back((double)obs[i]);
        obs_multi.push_back(v);
        obs_weights.push_back(1.0);
    }
    return this->fit(obs_multi, obs_weights);
}

double mixtureModel::fit(const vector<int>& obs, vector<double>& obs_weights){
    vector<vector<double> > obs_multi;
    for (int i = 0; i < obs.size(); ++i){
        vector<double> v;
        v.push_back((double)obs[i]);
        obs_multi.push_back(v);
    }
    return this->fit(obs_multi, obs_weights);
}

double mixtureModel::fit(const vector<vector<int> >& obs){
    vector<vector<double> > obsf;
    vector<double> obs_weights;
    for (int i = 0; i < obs.size(); ++i){
        vector<double> v;
        for (int j = 0; j < obs[i].size(); ++j){
            v.push_back((double)obs[i][j]);
        }
        obsf.push_back(v);
        obs_weights.push_back(1.0);
    }
    return this->fit(obsf, obs_weights);
}

double mixtureModel::fit(const vector<vector<int> >& obs, vector<double>& obs_weights){
    vector<vector<double> > obsf;
    for (int i = 0; i < obs.size(); ++i){
        vector<double> v;
        for (int j = 0; j < obs[i].size(); ++j){
            v.push_back((double)obs[i][j]);
        }
        obsf.push_back(v);
    }
    return this->fit(obsf, obs_weights);
}



double mixtureModel::fit(const vector<vector<double> >& obs){
    vector<double> obs_weights;
    for (int i = 0; i < obs.size(); ++i){
        obs_weights.push_back(1.0);
    }
    return this->fit(obs, obs_weights);
}

double mixtureModel::fit(const vector<vector<double> >& obs, vector<double>& obs_weights){
    
    if (this->n_components < 0){
        fprintf(stderr, "ERROR: model not initialized\n");
        exit(1);
    }
    if (obs.size() != obs_weights.size()){
        fprintf(stderr, "ERROR: fitting on %ld observations, but %ld observation weights provided.\n",
            obs.size(), obs_weights.size());
        exit(1);
    }
    int ni = this->dists[0].get_n_inputs();
    for (int i = 1; i < this->n_components; ++i){
        if (this->dists[i].get_n_inputs() != ni){
            fprintf(stderr, "ERROR: different number of inputs for component distributions:\n");
            fprintf(stderr, "dist 0: %d vs dist %d: %d\n", ni, i, this->dists[i].get_n_inputs());
            exit(1);
        }
    }

    // How many observations in data set?
    long int n_sites = obs.size();
    
    if (n_sites == 0){
        fprintf(stderr, "ERROR: cannot fit mixture model on 0 sites\n");
        exit(1);
    }
    // Is every component distribution frozen -- i.e. no need to update?
    bool all_frozen = true;
    for (int j = 0; j < this->n_components; ++j){
        if (e_only){
            this->dists[j].frozen = true;
        }
        if (!this->dists[j].frozen){
            all_frozen = false;
            break;
        }
    }
    
    // How many component distributions?
    if (n_components <= 1){
        return this->fit_single(obs, obs_weights);
    }
    // Need to create "responsibility matrix," storing the likelihoods of 
    // observations under each distribution multiplied by distribution weight
    // and normalized so that rows sum to 1. We make this accessible to outside
    // users, but need to remember to properly initialize and free it.
    this->init_responsibility_matrix(n_sites);
    
    // What are the stopping criteria?
    double delta_thresh = this->delta_thresh;
    int max_its = this->maxits;

    // Useful reference: https://stephens999.github.io/fiveMinuteStats/intro_to_em.html
    
    double delta = 99;
    int its = 0;
    double loglik_prev = 0.0;

    // Keep track of the mean of each dimension of each observation under each component dist
    double member_weight_sums[n_components];
    double** mean_sums = new double*[n_components];
    for (int j = 0; j < n_components; ++j){
        mean_sums[j] = new double[obs[0].size()]; 
    }
    
    // How small can the variance be before the component is zeroed out? 
    float var_thresh = 0.0;

    // Has the model reached a point where it only contains one component?
    bool one_component = false;
    
    // Make observation weights sum to num observations (this will keep weights from throwing off
    // BIC calculation)
    double obs_weight_sum = 0.0; 
    for (int i = 0; i < obs_weights.size(); ++i){
        obs_weight_sum += obs_weights[i];
    }
    // Try to prevent underflow by scaling things up by this number
    // (will divide later)
    double obs_weight_scale = (double)obs_weights.size();
    for (int i = 0; i < obs_weights.size(); ++i){
        //obs_weights[i] /= obs_weight_sum;
        if (obs_weights[i] > 0){
            obs_weights[i] = pow(2, log2(obs_weights[i]) - log2(obs_weight_sum) + log2(obs_weight_scale));
            if (isnan(obs_weights[i])){
                fprintf(stderr, "weight nan! %f %f %f\n", obs_weights[i], obs_weight_sum, obs_weight_scale);
                exit(1);
            }
        }
    }
    while (delta > delta_thresh && its < max_its && !one_component){
        
        // E - step: determine membership weights. 
        // Keep track of means of each component of each observation
        for (int j = 0; j < n_components; ++j){
            member_weight_sums[j] = 0.0;
            for (int k = 0; k < obs[0].size(); ++k){
                mean_sums[j][k] = 0.0;
            }
        }
        for (int i = 0; i < n_sites; ++i){
            if (obs_weights[i] == 0.0){
                continue;
            }
            double rowsum = 0;
            // Prevent underflow by dividing densities by the smallest density encountered
            double llmax = 0;
            double ll[n_components];
            for (int j = 0; j < n_components; ++j){
                if (this->weights[j] > 0.0){
                    ll[j] = log2(weights[j]) + this->dists[j].loglik(obs[i]);
                    if (isnan(ll[j]) || isinf(ll[j])){
                        fprintf(stderr, "Invalid log likelihood returned by distribution!\n");
                        fprintf(stderr, "weight %f\n", weights[j]);
                        fprintf(stderr, "Input:\n");
                        for (int k = 0; k < obs[i].size(); ++k){
                            fprintf(stderr, "%f ", obs[i][k]);
                        }
                        fprintf(stderr, "\n");
                        fprintf(stderr, "Distribution:\n");
                        this->dists[j].print(0);
                        exit(1);
                    }
                    if (llmax == 0 || ll[j] > llmax){
                        llmax = ll[j];
                    }
                }
                else{
                    ll[j] = 0.0;
                }
            }
            for (int j = 0; j < n_components; ++j){
                if (this->weights[j] > 0.0){
                    this->responsibility_matrix[i][j] = pow(2, ll[j] - llmax);
                    rowsum += this->responsibility_matrix[i][j];
                }
                else{
                    this->responsibility_matrix[i][j] = 0.0;
                }
            }
            if (rowsum == 0){
                rowsum = 1.0;
            }

            if (its == 0){
                // Compute log likelihood for initial fit.
                loglik_prev += log2(rowsum) + llmax;    
            }
            for (int j = 0; j < n_components; ++j){
                // Normalize responsibility matrix entries for this row
                this->responsibility_matrix[i][j] /= rowsum;

                // Track weight sums for each component
                member_weight_sums[j] += this->responsibility_matrix[i][j] * obs_weights[i];
                if (isnan(member_weight_sums[j])){
                    fprintf(stderr, "mws nan %d) j\n", j);
                    fprintf(stderr, "w %f rm %f rs %f\n", obs_weights[i], responsibility_matrix[i][j],
                        rowsum);
                    exit(1);
                } 

                for (int k = 0; k < obs[i].size(); ++k){
                    //mean_sums[j][k] += pow(2, log2(this->responsibility_matrix[i][j]) + 
                    //    log2(obs[i][k]) + log2(obs_weights[i]));
                    mean_sums[j][k] += this->responsibility_matrix[i][j] * obs[i][k] * obs_weights[i];
                }
            }
        }

        // M-step: update parameters
        for (int j = 0; j < n_components; ++j){
            if (weights[j] > 0){
                // Distribution weight
                // This is already normalized by including the (normalized) observation weights
                // in the sum
                double new_weight = member_weight_sums[j] / obs_weight_scale;
                if (isnan(new_weight)){
                    new_weight = 0.0;
                }
                this->weights[j] = new_weight;
            
                if (!e_only && new_weight > 0){
                    bool can_update = true;

                    // Compute summary statistics (mean and variance of each dimension of observations)
                    vector<double> means;
                    vector<double> vars;
                    map<pair<int, int>, double> covs;
                                        
                    for (int k = 0; k < obs[0].size(); ++k){
                        double mean = mean_sums[j][k] / member_weight_sums[j];
                        
                        if (dists[j].needs_cov){
                            // We must compute covariances with other observations
                            for (int l = k + 1; l < obs[0].size(); ++l){
                                pair<int, int> key = make_pair(k, l);
                                covs.insert(make_pair(key, 0.0));
                            }
                        }

                        double var = 0.0;
                        for (int i = 0; i < obs.size(); ++i){
                            
                            // If observation equals mean, this will contribute nothing
                            // to variance and cause a NaN inside the logarithm
                            if (obs[i][k] != mean && this->responsibility_matrix[i][j] != 0 &&
                                obs_weights[i] != 0){
                                double varlog = log2(this->responsibility_matrix[i][j]) + 
                                    2*log2(abs(obs[i][k] - mean)) + log2(obs_weights[i]);
                                var += pow(2, varlog);
                            }

                            if (dists[j].needs_cov){
                                for (int l = k + 1; l < obs[0].size(); ++l){
                                    pair<int, int> key = make_pair(k, l);
                                    double meanl = mean_sums[j][l] / member_weight_sums[j];
                                    if (obs[i][k] != mean && obs[i][l] != meanl && 
                                        this->responsibility_matrix[i][j] != 0 && 
                                        obs_weights[i] != 0){
                                        // Can't take logarithm since we might have negative numbers
                                        // here (not squaring)
                                        double covar = this->responsibility_matrix[i][j] * 
                                            (obs[i][k] - mean)*(obs[i][l] - meanl) * obs_weights[i];
                                        covs[key] += covar;
                                    }
                                }
                            }

                        }
                        var /= weights[j];
                        var /= obs_weight_scale;
                        
                        if (dists[j].needs_cov){
                            for (map<pair<int, int>, double>::iterator c = covs.begin(); 
                                c != covs.end(); ++c){
                                c->second /= weights[j];
                                c->second /= obs_weight_scale;
                            }
                        }

                        if (isnan(mean) || isinf(mean)){
                            fprintf(stderr, "mean %f %f\n", mean_sums[j][k], member_weight_sums[j]);
                            fprintf(stderr, "%f\n", weights[j]);
                            exit(1);
                        }
                        if (isnan(var) || isinf(var)){
                            fprintf(stderr, "var %f %f %f\n", var, weights[j], obs_weight_scale);
                            exit(1);
                        }
                        means.push_back(mean);
                        vars.push_back(var);
                        if (var <= 1e-8){
                            can_update = false;
                        }
                    }
                    
                    if (can_update){
                        // Allow distribution to update using method of moments
                        can_update = this->dists[j].update( means, vars, covs );
                        if (!can_update){
                            // Created illegal parameter values. Eliminate this
                            // distribution.
                            this->weights[j] = 0.0;
                            
                            // Give remaining weight to other components
                            double weightsum = 0.0;
                            for (int j2 = 0; j2 < n_components; ++j2){
                                weightsum += this->weights[j2];
                            }
                            for (int j2 = 0; j2 < n_components; ++j2){
                                this->weights[j2] /= weightsum;
                            }
                        }
                    }
                }
            } 
        }
        // Let external functions hook into the update process
        if (this->has_callback_fun){
            this->callback_fun(*this, this->shared_params);
            /*
            vector<mixtureDist*> distvec;
            for (int i = 0; i < this->dists.size(); ++i){
                distvec.push_back(&this->dists[i]);
            }
            this->callback_fun(distvec, this->weights, this->shared_params);
            */
        }
        
        for (int j = 0; j < this->weights.size(); ++j){
            if (this->weights[j] == 1.0){
                one_component = true;
                break;
            }
        }
        double loglik = this->compute_loglik(obs, obs_weights); 
        if (isinf(loglik)){
            fprintf(stderr, "LL inf\n");
            exit(1);
        }
        if (this->print_lls){
            fprintf(stderr, "LL %.3f -> %.3f (Improvement: %.3f)\n", loglik_prev/log2(exp(1)), 
                loglik/log2(exp(1)), (loglik-loglik_prev)/log2(exp(1)));
            if (this->print_dists){
                fprintf(stderr, "\n");
                this->print();
                fprintf(stderr, "\n");
            }
        }

        //delta = abs(loglik-loglik_prev);
        if (all_frozen){
            delta = 0.0;
        }
        else{
            delta = loglik-loglik_prev;
            if (isnan(delta)){
                fprintf(stderr, "delta NAN %f %f\n", loglik_prev, loglik);
                this->print();
                exit(1);
            }
        }

        loglik_prev = loglik;
        ++its;

    }
    
    if (this->print_lls){
        if (delta <= delta_thresh){
            fprintf(stderr, "Converged in %d iterations\n", its);
        }
        else{
            fprintf(stderr, "Did not converge!\n");
            fprintf(stderr, "\titeration %d max %d\n", its, max_its);
            fprintf(stderr, "\tdelta %.3f thresh %.3f\n", delta, delta_thresh);
            fprintf(stderr, "\tweights:\n");
            for (int j = 0; j < this->weights.size(); ++j){
                fprintf(stderr, "\t\t%.3f\n", this->weights[j]);
            }
        }
    }
    // Assign observations to likeliest component of origin.
    for (int i = 0; i < n_sites; ++i){
        double maxweight = 0;
        int max_component = -1;
        for (short j = 0; j < n_components; ++j){
            if (this->responsibility_matrix[i][j] > maxweight){
                maxweight = this->responsibility_matrix[i][j];
                max_component = j;
            }
        }
        this->assignments.push_back(max_component);
    }
    for (int j = 0; j < n_components; ++j){
        delete[] mean_sums[j];
        mean_sums[j] = NULL;
    }
    delete[] mean_sums;
    // Store log likelihood of fit model
    this->loglik = loglik_prev;
    // Compute BIC and AIC
    this->compute_bic(obs.size() * obs[0].size());    
    return loglik_prev;
   
}

void mixtureModel::trigger_callback(){
    if (this->has_callback_fun){
        this->callback_fun(*this, this->shared_params);
        /*
        vector<mixtureDist*> distvec;
        for (int i = 0; i < this->dists.size(); ++i){
            distvec.push_back(&this->dists[i]);
        }
        this->callback_fun(distvec, this->weights, this->shared_params);
        */
    }
}

double mixtureModel::compute_loglik(const vector<vector<double> >& obs, vector<double>& obs_weights){
    double loglik = 0;

    int n_components_pos = 0;
    for (int j = 0; j < this->n_components; ++j){
        if (this->weights[j] > 0.0){
            n_components_pos++;
            if (n_components_pos > 1){
                break;
            }
        }
    }
    
    bool one_component = n_components_pos == 1;
    
    for (int i = 0; i < obs.size(); ++i){
        if (obs_weights[i] == 0.0){
            continue;
        }
        double ll_row[n_components];
        double ll_row_max = 0.0;
        for (int j = 0; j < n_components; ++j){
            ll_row[j] = 0;
            if (this->weights[j] > 0){
                if (!one_component){
                    ll_row[j] += log2(this->weights[j]);
                }
                ll_row[j] += this->dists[j].loglik( obs[i] );
                if (ll_row_max == 0.0 || ll_row[j] > ll_row_max){
                    ll_row_max = ll_row[j];
                }
            }
        }
        double row_sum = 0.0;
        for (int j = 0; j < this->n_components; ++j){
            if (this->weights[j] > 0){
                row_sum += pow(2, ll_row[j] - ll_row_max);
            }
        }
        loglik += log2(row_sum) + ll_row_max;
        if (isinf(loglik)){
            fprintf(stderr, "LL inf\n");
            fprintf(stderr, "row sum %f\n", row_sum);
            fprintf(stderr, "lrowsum %f\n", log2(row_sum));
            fprintf(stderr, "max %f\n", ll_row_max);
            this->print();
            exit(1);
        }
    }
    return loglik;
}
void mixtureModel::compute_bic(long int n_obs){
    if (this->loglik == 0.0){
        fprintf(stderr, "ERROR: log likelihood not set. Fit model before computing BIC/AIC.\n");
        exit(1);
    }

    // Get number of parameters
    int n_params = 0;
    for (int i = 0; i < this->n_components; ++i){
        // How many parameters does this distribution have? 
        n_params += this->dists[i].get_n_parameters();
    }

    // Each distribution's weight is also a parameter
    if (this->n_components > 1){
        n_params += (int)this->weights.size();
    }
    
    this->bic = (float)n_params * log(n_obs) - 2.0 * (this->loglik / log2(exp(1)));
    this->aic = 2.0 * (float)n_params - 2.0 * (this->loglik / log2(exp(1)));
}

double mixtureModel::fit_single(const vector<vector<double> >& obs, vector<double>& obs_weights){
    
    if (this->dists.size() == 0){
        fprintf(stderr, "ERROR: to fit model, please initialize parameters for a single model \
with dummy parameters. This is necessary to judge how many parameters are being fit (values are irrelevant.\n");
        exit(1);
    }
    
    // In this case, the responsibility matrix will just be a column of 1s
    this->init_responsibility_matrix(obs.size());
    for (int i = 0; i < obs.size(); ++i){
        this->responsibility_matrix[i][0] = 1.0;
    }

    // Ignore component weights; find best fit single distribution. Fill in log likelihood,
    // AIC, and BIC.
    
    vector<double> means;
    vector<double> vars;
    map<pair<int, int>, double> covs;

    double obs_weight_sum = 0.0;
    for (int i = 0; i < obs_weights.size(); ++i){
        obs_weight_sum += obs_weights[i];
    }
    double obs_weight_scale = (double)obs_weights.size();
    for (int i = 0; i < obs_weights.size(); ++i){
        if (obs_weights[i] > 0.0){
            obs_weights[i] = pow(2, log2(obs_weights[i]) - log2(obs_weight_sum) + log2(obs_weight_scale));
        }
    }

    for (int j = 0; j < obs[0].size(); ++j){
        means.push_back(0.0);
        vars.push_back(0.0);
        if (dists[0].needs_cov){
            for (int k = j + 1; k < obs[0].size(); ++k){
                pair<int, int> key;
                covs.insert(make_pair(key, 0.0));
            }
        }
    }
    double denom = (double)obs.size();
    for (int i = 0; i < obs.size(); ++i){
        if (obs_weights[i] == 0.0){
            continue;
        }
        for (int j = 0; j < obs[i].size(); ++j){
            means[j] += obs_weights[i] * obs[i][j];
        }
    }
    
    for (int j = 0; j < means.size(); ++j){
        means[j] /= obs_weight_scale;
    }
    
    for (int i = 0; i < obs.size(); ++i){
        if (obs_weights[i] == 0.0){
            continue;
        }
        for (int j = 0; j < obs[i].size(); ++j){
            vars[j] += pow(obs[i][j] - means[j], 2) * obs_weights[i];
            if (dists[0].needs_cov){
                for (int k = j + 1; k < obs[i].size(); ++k){
                    pair<int, int> key = make_pair(j, k);
                    covs[key] += (obs[i][j] - means[j])*(obs[i][k] - means[k]) * obs_weights[i];
                }
            }
        }
    }

    for (int j = 0; j < vars.size(); ++j){
        vars[j] /= obs_weight_scale;
    }
    for (map<pair<int, int>, double>::iterator cov = covs.begin(); cov != covs.end(); ++cov){
        cov->second /= obs_weight_scale;
    }

    // Fit the model (only one distribution)
    this->dists[0].update( means, vars, covs );
    
    this->weights.clear();
    this->weights.push_back(1.0);

    // Compute the log likelihood of the model
    double loglik = compute_loglik(obs, obs_weights);
    this->loglik = loglik;
    if (loglik == 0.0){
        fprintf(stderr, "LL == 0?\n");
        for (int i = 0; i < obs_weights.size(); ++i){
            fprintf(stderr, "w %d) %f\n", i, obs_weights[i]);
        }
        exit(1);
    }
    this->compute_bic(obs.size() * obs[0].size());
    return loglik;
}


void mixtureModel::print() const{
    fprintf(stderr, "=== Mixture model with %d components ===\n", this->n_components);
    bool needs_newline = false;
    if (this->loglik != 0){
        fprintf(stderr, "LL = %.3f ", this->loglik);
        needs_newline = true;
    }
    if (this->bic != 0 || this->aic != 0){
        fprintf(stderr, "BIC %.3f AIC %.3f", this->bic, this->aic);
        needs_newline = true;
    }
    if (needs_newline){
        fprintf(stderr, "\n");
    }
    //fprintf(stderr, "log2 likelihood = %.3f BIC %.3f AIC %.3f\n", this->loglik,
    //    this->bic, this->aic);
    for (int i = 0; i < this->n_components; ++i){
        fprintf(stderr, "   Component weight: %.3f\n", this->weights[i]);
        this->dists[i].print(2);
    }
}

