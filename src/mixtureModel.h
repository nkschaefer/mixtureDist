#ifndef _MIXTUREMODEL_H
#define _MIXTUREMODEL_H
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

// MIXTURE MODELS =====

// Generalized function to handle globally shared parameters of component distributions.
// Runs just after all distributions have updated themselves (by fitting to the data).
// Gives pointers to all distributions in a set with shared parameters, and those 
// distributions' weights (scaled so that they sum to one within this set), so that an
// external program can reconcile these parameters based on some other, global 
// parameters they might depend on. It is up to that program to decide what to do and
// how, or if, to adjust each distribution's parameters.

class mixtureModel;

typedef std::function< void( mixtureModel&,
    std::vector<double>& ) > callback;

//typedef std::function< void( const std::vector<mixtureDist*>&,
//    const std::vector<double>&,
//    std::vector<double>& ) > callback;

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
        std::vector<double> shared_params;

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
        
        // Allow users to hook into update after the M step with an external function
        bool set_callback(callback fun);
        bool set_callback(callback fun, std::vector<double> shared_meta_params);
        void trigger_callback();

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
        // What function should be called (externally) to reconcile distributions that share
        // parameters, after each maximization step?
        callback callback_fun;
        bool has_callback_fun;
        
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
