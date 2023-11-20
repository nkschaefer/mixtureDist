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
#include "common.h"
#include "incbeta.h"
#include "cdflib.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Uses Stirling's approximation
 * returns log2 factorial
 */
double factorial(int n){
    if (n == 0 || n == 1){
        return 0;      
    }
    else{
       return 0.5*(log2(2) + log2(M_PI) + log2(n)) + (float)n*(log2(n) - log2(exp(1)));
    }
}

/**
 *  Log2 PMF of Poisson distribution at x with mean l
 */
double pois2(double x, double l){
    if (l < 1){
        return -1;
    }
    double xfac = factorial(x);
    return x * log2(l) + -((float)l/log(2)) - xfac;
}

/**
 *  Log2 PMF of Poisson distribution at x with mean l
 */
double pois(int x, int l){
    if (l < 1){
        return -1;
    }
    double xfac = factorial(x);
    /*
    double xfac;
    if (x == 0){
        xfac = 0;    
    }
    else{
        static map<int, float> xfac_cache;
        if (xfac_cache.count(x) > 0){
            xfac = xfac_cache[x];
        }
        else{
            xfac = log2(x);
            int i = x-1;
            while (i >= 1){
                xfac = xfac + log2(i);
                i--;
            }
            xfac_cache.insert(make_pair(x, xfac));
        }
    }
    */
    return x * log2(l) + -((float)l/log(2)) - xfac;
}

double pois_cdf(int x, int l){
    
    double p;
    double q;
    double s = (double)x;
    double xlam = (double)l;
    cumpoi(&s, &xlam, &p, &q);
    return log2(p);
    /*
    int which = 1;
    double p = 0.0;
    double q = 0.0;
    double s = (double) x;
    double xlam = (double) l;
    int status;
    double bound;
    // which = which argument is to be calculated
    //  1: calculate p and q from s and xlam
    //  2: calculate a from p, q, and xlam
    //  3: calculate xlam from p, q, and s
    //  p = the cumulation from 0 to s of the poisson density (prob)
    //  q = 1-p
    //  s = upper limit of cumulation of poisson CDF
    //  xlam = mean of poisson distribution
    //  status:
    //   0: success
    //   -I: I parameter out of range
    //   +1: answer lower than lowest search bound
    //   +2: answer higher than greatest search bound
    //   +3: p + q != 1
    //   bound is defined if status nonzero
    //    status < 0: value exceeded by parameter i
    //    status 1 or 2: search bound that was exceeded.
    
    cdfpoi(&which, &p, &q, &s, &xlam, &status, &bound);
    
    if (status < 0){
        fprintf(stderr, "err computing poisson CDF\n");
    }   
    return log2(p) - log2(q);
    */
}

double dlgamma(double x, double a, double b){
    double term1 = a * log2(b);
    double term2 = lgammaf(a) / log(2);
    double term3 = (a - 1.0) * log2(x);
    double term4 = (-b * x) * log2(exp(1));
    return (term1-term2) + term3 + term4;
}

pair<double, double> gamma_moments(double mu, double var){
    double alpha = pow(mu, 2) / var;
    //double beta = var / mu;
    double beta = mu / var;
    return make_pair(alpha, beta);
}

/**
 * CDF of gamma distribution
 */
double pgamma(double x, double a, double b){
    // Convert from shape/rate to shape/scale parameterization
    double shape = a;
    //double scale = 1.0/b;
    double scale = b; // I think they confused "scale" for "rate" in documentation
    int status;
    int which = 1; // calculate p & q from x, shape, and scale
    double p;
    double q;
    double bound;
    cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
    if (status != 0){
        // -I = input parameter index I is out of range
        //  bound set to value exceeded by parameter I
        // 1 = answer lower than lowest search bound
        //  bound set to bound that was exceeded
        // 2 = answer higher than highest search bound
        //  bound set to bound that was exceeded
        // 3 p + q != 1
        // 10 gamma or inverse gamma routine didn't converge
        return -1.0;
    }
    return p;
}

/**
 * Log2 PDF of exponential distribution at x with mean l
 */
double expPDF(float x, float l){
    if (l <= 0 || x <= 0){
        return -1;
    }
    return log2(l) - l*x*log2(exp(1));
}

/**
 * Binomial coefficient
 */
int binom_coef(int n, int k){
    if (k > n){
        return 0;
    }
    if (n == 0){
        return 0;
    }
    if (k == 0){
        return 1;
    }
    if (n == k){
        return 1;
    }
    static map<int, float> fac_cache;
    
    float nfac = 0.0;
    if (n > 0){
        if (fac_cache.count(n) > 0){
            nfac = fac_cache[n];
        }
        else{
            for (int i = n; i >= 1; --i){
                nfac += log2(i);
            }
            fac_cache.insert(make_pair(n, nfac));
        }
    }
    float kfac = 0.0;
    if (k > 0){
        if (fac_cache.count(k) > 0){
            kfac = fac_cache[k];
        }
        else{
            for (int i = k; i >= 1; --i){
                kfac += log2(i);
            }
            fac_cache.insert(make_pair(k, kfac));
        }
    }
    float nminuskfac = 0.0;
    if (n - k > 0){
        if (fac_cache.count(n-k) > 0){
            nminuskfac = fac_cache[n-k];
        }
        else{
            for (int i = n - k; i >= 1; --i){
                nminuskfac += log2(i);
            }
            fac_cache.insert(make_pair(n-k, nminuskfac));
        }
    }
    float nchoosek = nfac - (kfac + nminuskfac);
    return pow(2.0, nchoosek);
    /*
    int nfac = n;
    for (int x = n - 1; x > 1; x--){
        nfac *= x;
    }
    int kfac = k;
    for (int x = k - 1; x > 1; x--){
        kfac *= x;
    }
    int nminuskfac = n - k;
    for (int x = n - k - 1; x > 1; x--){
        nminuskfac *= x;
    }
    return nfac / (kfac * nminuskfac);
    */
}
double binom_coef_log(double n, double k){
    if (k > n){
        return 0;
    }
    if (n == 0){
        return 0;
    }
    if (n == k){
        return 0;
    }
    if (k == 0){
        return 0; // 1 way to choose
    }
    double nf = (double)n;
    double kf = (double)k;
    
    /* 
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.binom.html
    float term1 = lgammaf(nf + 1.0);
    float term2 = lgammaf(kf + 1.0);
    float term3 = lgammaf(nf-kf + 1.0);
    float res = term1 / (term2 + term3);
    return res/M_LN2; 
    */

    // https://math.stackexchange.com/questions/64716/approximating-the-logarithm-of-the-binomial-coefficient
    // Use Stirling's approximation
    //return nf*log2(nf) - kf*log2(kf) - (nf-kf)*log2(nf-kf);

    return nf*log2(nf) - kf*log2(kf) - (nf-kf)*log2(nf-kf) + 0.5*(log2(nf) - log2(kf) - log2(nf - kf) - log2(2*M_PI));

    
    /*
    static map<int, float> fac_cache;
    
    float nfac = 0.0;
    if (n > 0){
        if (fac_cache.count(n) > 0){
            nfac = fac_cache[n];
        }
        else{
            for (int i = n; i >= 1; --i){
                nfac += log2(i);
            }
            fac_cache.insert(make_pair(n, nfac));
        }
    }
    float kfac = 0.0;
    if (k > 0){
        if (fac_cache.count(k) > 0){
            kfac = fac_cache[k];
        }
        else{
            for (int i = k; i >= 1; --i){
                kfac += log2(i);
            }
            fac_cache.insert(make_pair(k, kfac));
        }
    }
    float nminuskfac = 0.0;
    if (n - k > 0){
        if (fac_cache.count(n-k) > 0){
            nminuskfac = fac_cache[n-k];
        }
        else{
            for (int i = n - k; i >= 1; --i){
                nminuskfac += log2(i);
            }
            fac_cache.insert(make_pair(n-k, nminuskfac));
        }
    }
    float nchoosek = nfac - (kfac + nminuskfac);
    return nchoosek;
    */
}

float fac_log(int n){
    float nfac = 0.0;
    static map<int, float> fac_cache;
    if (n > 0){
        if (fac_cache.count(n) > 0){
            nfac = fac_cache[n];
        }
        else{
            for (int i = n; i >= 1; --i){
                nfac += log2(i);
            }
            fac_cache.insert(make_pair(n, nfac));
        }
    }
    return nfac;
}

float fac_gamma(int n){
    return fac_log(n-1);
}

/**
 * Only valid for integers!
 *
float betafunc(int a, int b){
    float term1 = log2(a + b) - log2(a*b);
    float term2 = binom_coef_log(a + b, a);
    return term1 - term2;
}
*/

/**
 * Log2 of the beta function
 * Valid for all real numbers
 */
float lbetaf(float a, float b){
    return ((lgammaf(a) + lgammaf(b)) - lgammaf(a + b))/M_LN2;
}

/**
 * Regular (non-log) beta function
 */
long double betaf(float a, float b){
    long double term1 = tgammal(a) / tgammal(a+b);
    long double term2 = tgammal(b) / tgammal(a+b);
    return term1*term2;
}

/**
 * log2 transformed
 * x = number successes
 * n = number trials
 * a = alpha for beta dist
 * b = beta for beta dist
 */
double dbetabin(double x, double n, double a, double b){
    double term1 = binom_coef_log(n, x);
    double term2 = lbetaf(x + a, (n - x) + b);
    double term3 = lbetaf(a, b);
    return term1 + term2 - term3;
}

/**
 * Log PDF of beta distribution
 */
float dlbeta(float x, float a, float b){
    float term1 = (a-1) * log2(x);
    float term2 = (b-1) * log2(1.0-x);
    //float term3 = betafunc(a, b);
    float term3 = lbetaf(a, b);
    return (term1 + term2) - term3;
}

/**
 * non-log PDF of beta distribution
 */
float dbeta(float x, float a, float b){
    if (x == 0.0 || x == 1.0){
        return 0.0;
    }
    else{
        return pow(2, dlbeta(x, a, b));
    }
}

/** 
 * CDF of beta distribution
 */
float pbeta(float x, float a, float b){
    return incbeta(a, b, x); 
}

/**
 * Fit a beta distribution by method of moments
 */
pair<float, float> betafit(vector<float>& vals){
    pair<float, float> muvar = welford(vals);
    float mean = muvar.first;
    float var = muvar.second;
    float alpha = mean*((mean*(1-mean))/var - 1);
    float beta = (alpha * (1-mean))/mean;
    return make_pair(alpha, beta);
}

/**
 * Finds best fit scale parameter for a beta distribution with fixed
 * location parameter
 */
pair<float, float> betafit_disp(vector<float>& vals, float loc){
    float varsum = 0.0;
    for (int i = 0; i < vals.size(); ++i){
        varsum += pow(vals[i]-loc, 2);
    }
    varsum /= ((float)vals.size()-1);
    float alpha = loc*((loc*(1-loc))/varsum -1);
    float beta = (alpha * (1-loc))/loc;
    return make_pair(alpha, beta);
}

float dnbinom_mup(int x, float mu, float p){
    // mu is mean
    // p is mean/variance
    float rf = pow(mu,2)/((mu/p) - mu);
    int r = (int)round(rf);

    float term1 = binom_coef_log(x + r - 1, x);
    float term2 = rf * log2(p);
    float term3 = (float)x * log2(1.0 - p);
    
    fprintf(stderr, "x %d mu %f p %f | %f %f %f\n", x, mu, p, term1, term2, term3);
    fprintf(stderr, "rf %f r %d\n", rf, r);
    fprintf(stderr, "%f %f | %f\n", pow(mu, 2), mu/p, mu/p - mu);
    return term1 + term2 + term3;
    /*
    float rf = mu - log2(var - mu);
    int r = (int)round(pow(rf, 2));
    //int r = (int)round(pow(mu, 2) / (var - mu));
    float p = mu/var;
    float term1 = binom_coef_log(x + r - 1, r - 1); 
    float term2 = x * log2(1.0-p);
    float term3 = r * log2(p);
    if (isnan(term1) || isnan(term2) || isnan(term3)){
        fprintf(stderr, "%d %f %f\n", x, mu, var);
        fprintf(stderr, "%f %f %f\n", term1, term2, term3);
        exit(1);
    }
    return term1 + term2 + term3;
    */
}

/**
 * Log negative binomial PDF of r successes in k trials with success probability p
 *
 * phi = mean^2/(variance - mean)
 *  var*phi - mu*phi = mu^2
 *  var*phi = mu(mu + phi)
 *  var = mu(mu+phi)/phi
 *  mean/variance = mu/(mu(mu+phi)/phi) = (mu*phi)/(mu*(mu+phi))
 *      = phi/(mu+phi)
 *  
 * To convert between these parameters and the Wikipedia definition:
 * x = k (user input)
 * p = mean/variance = phi/(mu+phi)
 * r = phi
 *
 */
float dnbinom(int x, int mu, float phi){
    //return binom_coef(k + r - 1, k) + (k * log2(1.0-p)) + r * log2(p);
    float term1 = binom_coef_log(x + phi - 1, x);
    float term2 = x * (log2(mu) - log2(mu + phi));
    float term3 = phi * (log2(phi) - log2(mu + phi));
    return term1 + term2 + term3;
}

/**
 * Negative binomial CDF 
 * uses mu/phi parameterization
 */
float pnbinom(int x, int mu, float phi){
    float p = phi/((float)mu+phi);
    // incbeta syntax is a, b, x for I_x(a, b)
    return incbeta(phi, x+1, p); 
}

/**
 * Fit negative binomial distribution using method of moments
 */
pair<int, float> nbinom_moments(float mean, float var){
    int mu = (int)round(mean);
    float phi = pow((float)mean, 2) / (var - mean);
    return make_pair(mu, phi);
}

/**
 * Log binomial PDF of k successes in n trials with success probability p
 */
double dbinom(double n, double k, double p){
    double nchoosek = binom_coef_log(n, k);
    double res =  nchoosek + k*log2(p) + (n-k)*log2(1-p);
    
    return res;
}

/**
 * Log multinomial PDF
 */
double dlmultinom(const vector<double>& x, const vector<double>& p){
    if (x.size() != p.size()){
        fprintf(stderr, "ERROR: dlmultinom: length of x != length of p\n");
        exit(1);
    }
    double xsum = 1;
    double term2 = 0;
    double term3 = 0;
    for (int i = 0; i < x.size(); ++i){
        xsum += x[i];
        term2 += lgammaf(x[i] + 1);
        term3 += x[i] * log(p[i]);
    }
    double term1 = lgammaf(xsum);
    return (term1 - term2 + term3)/log(2);
}

/**
 * Log Dirichlet PDF
 */
double dldirichlet(const vector<double>& x, const vector<double>& alpha){
    if (x.size() != alpha.size()){
        fprintf(stderr, "ERROR: dldirichlet: length of x != length of alpha\n");
        exit(1);
    }
    double alpha_0 = 0.0;
    double xsum = 0.0;
    for (int i = 0; i < alpha.size(); ++i){
        alpha_0 += alpha[i];
        xsum += x[i];
    }
    if (xsum != 1.0){
        fprintf(stderr, "ERROR: dldirichlet: x vector does not sum to 1 (%f)\n", xsum);
        exit(1);
    }
    double betaprod = 0.0;
    double betadenom = lgammaf(alpha_0);
    double term2 = 0.0;
    for (int i = 0; i < alpha.size(); ++i){
        betaprod += lgammaf(alpha[i]);
        term2 += (alpha[i]-1) * log(x[i]);
    }
    double beta = betaprod - betadenom;
    return (term2 - beta)/log(2);
}

/** 
 * PMF of hypergeometric distribution
 * N population size
 * K success states in population
 * n sample size
 * k sample successes
 */
double pmfhyper(int N, int K, int n, int k){
    return (binom_coef_log(K, k) + binom_coef_log(N-K, n-k)) - binom_coef_log(N,n);    
}
double binom_cdf_lower(int n, int k, float p){
    // Regularized incomplete beta function
    // = B(x,a,b)
    // x = q = 1 - p
    // a = n-k
    // b = 1 + k
    return incbeta(n-k, 1+k,1-p);
}

/**
 * Log PDF of normal distribution
 */
float dnorm(float x, float mu, float sigma){
    static const float sqrt_2pi = 1.325748;
    float a = (x - mu) / sigma;
    
    return (-0.5*pow(a,2)) / log(2) - log2(sigma) - log2(sqrt(2*3.14159265358979));
    //return - log2(sigma) - sqrt_2pi - pow(a,2)/log(2);
}
float cnorm(float x, float mu, float sigma){
    //return log2(0.5) + log2(erfc(-((x-mu)/sigma) * M_SQRT1_2));
    return 0.5 * erfc(-((x-mu)/sigma) * M_SQRT1_2);
}

/**
 * Normal (log) approximation to the binomial distribution
 * (faster to compute)
 */
float binom_approx(int n, int k, float frac){
    float mu = (float)n * frac;
    float sigma = sqrt((float)n * frac * (1.0-frac));
    return dnorm(k, mu, sigma);
}

float binom_approx_scaled(int n, int k, float scale, float frac){
    float mu = (float)n * frac;
    float sigma = sqrt(scale*(float)n * frac * (1.0-frac));
    return dnorm(k, mu, sigma);
}

/**
 * Non-log PDF of normal distribution with params mu and sigma
 */
float gaussPDF(float x, float mu, float sigma){
    //return ((1.0/(sigma*sqrt(2*3.14159265358979)))*exp(-0.5*(pow((x-mu)/sigma, 2))));
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - mu) / sigma;
    return inv_sqrt_2pi / sigma * exp(-0.5f * a * a);
}

/**
 * Finds the point of intersection between two gaussians.
 */
float int2gauss(float mu1, float mu2, float sigma1, float sigma2){
    if (mu2 < mu1){
        // Flip.
        float tmp = mu2;
        float tmp2 = sigma2;
        mu2 = mu1;
        sigma2 = sigma1;
        mu1 = tmp;
        sigma1 = tmp2;
    }
    float sig21 = sigma1*sigma1;
    float sig22 = sigma2*sigma2;
    float mu21 = mu1*mu1;
    float mu22 = mu2*mu2;
    
    float a = (-1.0/sig21 + 1.0/sig22)/2;
    float b = -mu2/sig22 + mu1/sig21;
    float c = (mu22/sig22 - mu21/sig21 + log(sig22/sig21))/2;
    
    // Solve ax^2 + bx + c
    float sol1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
    float sol2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
    if (sol1 < 0 && sol2 > 0){
        return sol2;
    }
    else if (sol2 < 0 && sol1 > 0){
        return sol1;
    }
    else if (sol1 < sol2){
        return sol1;
    }
    return sol2;
    
}

/**
 * Expectation maximization fitting of normal mixture model to data.
 */
void gaussEM(vector<float>& obs, 
    vector<float>& weights, 
    vector<float>& mu, 
    vector<float>& sigma,
    vector<short>& assignments){

    if (weights.size() != mu.size() || weights.size() != sigma.size()){
        fprintf(stderr, "ERROR: weights, mu, and sigma must be the same size.\n");
        return;
    }
    /* 
    for (int i = 0; i < obs.size(); ++i){
        fprintf(stderr, "OBS %f\n", obs[i]);
        if (isnan(obs[i])){
            exit(1);
        }
    }
    */
    int n_components = mu.size();
    int n_sites = obs.size();
    float member_weights[n_sites][n_components];

    float delta_thresh = 0.001;
    int max_its = 1000;
    
    // Useful reference: https://stephens999.github.io/fiveMinuteStats/intro_to_em.html

    float delta = 99;
    int its = 0;
    while (delta > delta_thresh && its < max_its){
        
        // E - step: determine membership weights.
        
        float member_weight_sums[n_components];
        float mu_sums[n_components];
        //fprintf(stderr, "COMPONENTS ");
        for (int j = 0; j < n_components; ++j){
            member_weight_sums[j] = 0.0;
            mu_sums[j] = 0.0;
            //fprintf(stderr, "w %f m %f s %f ", weights[j], mu[j], sigma[j]);
        }
        //fprintf(stderr, "\n");
        for (int i = 0; i < n_sites; ++i){
            float rowsum = 0;
            for (int j = 0; j < n_components; ++j){
                if (weights[j] > 0.0){
                    float gpdf;
                    if (sigma[j] > 0){
                        gpdf = gaussPDF(obs[i], mu[j], sigma[j]);
                    }
                    else{
                        gpdf = 0.0;
                    }
                    member_weights[i][j] = weights[j] * gpdf;
                    rowsum += member_weights[i][j];
                }
                else{
                    member_weights[i][j] = 0.0;
                }
            }
            if (rowsum == 0){
                rowsum = 1.0;
            }
            for (int j = 0; j < n_components; ++j){
                
                member_weights[i][j] /= rowsum;
                //fprintf(stderr, "%f ", member_weights[i][j]);
                member_weight_sums[j] += member_weights[i][j];
                mu_sums[j] += member_weights[i][j] * obs[i];
            }
            //fprintf(stderr, "\n");
        }
        
        float maxdelta = -1;
        
        // M-step: update parameters
        for (int j = 0; j < n_components; ++j){
            if (weights[j] > 0){ 
                // Distribution weight
                float new_weight = member_weight_sums[j] / (float)n_sites;
                float delta_weight = abs(new_weight - weights[j]) / weights[j];
                if (delta_weight > maxdelta){
                    maxdelta = delta_weight;
                }
                
                weights[j] = new_weight;
                // Distribution mean
                float new_mu = mu_sums[j] / member_weight_sums[j];
                float delta_mu = abs(new_mu - mu[j]) / mu[j];
                if (delta_mu > maxdelta){
                    maxdelta = delta_mu;
                }
                mu[j] = new_mu;

                // Distribution sd
                float sigmasum = 0.0;
                for (int i = 0; i < obs.size(); ++i){
                    sigmasum += member_weights[i][j]*pow(obs[i] - new_mu, 2);
                }
                float new_sigma = sqrt(sigmasum / member_weight_sums[j]);
                //if (new_sigma == 0.0){
                //    new_sigma = 1e-4;
                //}
                float delta_sigma = abs(new_sigma - sigma[j]) / sigma[j];
                if (delta_sigma > maxdelta){
                    maxdelta = delta_sigma;
                }
                sigma[j] = new_sigma;
            } 
        }
           
        delta = maxdelta;

        ++its;
    }
    
    if (delta <= delta_thresh){
        fprintf(stderr, "Converged in %d iterations\n", its);
    }
    else{
        fprintf(stderr, "Reached maximum of %d iterations\n", its);
    }

    // At this step, weights, mu, and sigma have info for distributions
    // member_weights rows are sites - columns are probability site belongs to component j

    fprintf(stderr, "MIXTURE COMPONENTS\n");
    for (int j = 0; j < n_components; ++j){
        fprintf(stderr, "weight %f mu %f sd %f\n", weights[j], mu[j], sigma[j]);
        if (isnan(weights[j])){
            exit(1);
        }
        //if (isnan(weights[j]) || isnan(mu[j]) || isnan(sigma[j])){
        //    exit(1);
        //}
    }
    
    // Assign observations to likeliest component of origin.
    for (int i = 0; i < n_sites; ++i){
        float maxweight = 0;
        int max_component = -1;
        for (int j = 0; j < n_components; ++j){
            if (member_weights[i][j] > maxweight){
                maxweight = member_weights[i][j];
                max_component = j;
            }
        }
        assignments.push_back(max_component);
    }
}

/**
 * Convert location and scale parameters for beta distribution into
 * alpha and beta parameterization.
 */
pair<float, float> beta_loc_scale_to_a_b(float loc, float scale){
    //float alpha = (loc - scale*loc)/scale;
    //float beta = (1.0/scale) * (1-loc) + loc - 1;
    float alpha = loc*scale;
    float beta = scale - loc*scale;
    return make_pair(alpha, beta);
}

/**
 * Welford's algorithm for computing mean and variance in a single pass
 * Returns pair of (mean, variance)
 */
pair<double, double> welford(vector<float>& vals){
    double count = 0;
    double mean = 0;
    double M2 = 0;
    
    for (int i = 0; i < vals.size(); ++i){
        count++;
        double delta1 = vals[i]-mean;
        mean += delta1/count;
        double delta2 = vals[i]-mean;
        M2 += delta1*delta2;
    }
    return make_pair(mean, M2/(count-1.0));
}

// ===== mixtureDist =====

/**
 * Define log likelihood and update functions for pre-set probability distributions.
 */
double mixtureDist::ll_poisson(const vector<double>& input, 
    int start_idx,
    int n_inputs, 
    const vector<double>& params){
    return pois2(input[start_idx], params[0]);
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
    return dlgamma(input[start_idx], params[0], params[1]); 
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
    return dlmultinom(inputs_this, params);
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
}

mixtureDist::mixtureDist(string name, vector<double> params_input){
    this->n_components = 0;
    mixtureDist::auto_register();
    if (!this->add_dist(name, params_input)){
        exit(1);
    }
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
}

mixtureDist::mixtureDist(vector<string> names, vector<vector<double> > params_input){
    this->n_components = 0;
    mixtureDist::auto_register();
    for (int i = 0; i < names.size(); ++i){
        if (!this->add_dist(names[i], params_input[i])){
            exit(1);
        }
    }
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
}

mixtureDist::mixtureDist(vector<string> names, vector<double> weights){
    this->n_components = 0;
    mixtureDist::auto_register();
    for (int i = 0; i < names.size(); ++i){
        if (!this->add_dist(names[i], weights[i])){
            exit(1);
        }
    }
}

mixtureDist::mixtureDist(vector<string> names, vector<vector<double> > params_input, vector<double> weights){
    this->n_components = 0;
    mixtureDist::auto_register();
    for (int i = 0; i < names.size(); ++i){
        if (!this->add_dist(names[i], params_input[i], weights[i])){
            exit(1);
        }
    }
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
        fprintf(stderr, "%s%sDist with %d components:\n", indent_str.c_str(), namestr.c_str(),
            this->n_components);
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
}

void mixtureModel::init_responsibility_matrix(int n_obs){
    if (this->n_obs != -1){
        // Responsibility matrix was formerly initialized for a different
        // data set
        this->free_responsibility_matrix();
    }
    
    this->n_obs = n_obs;

    // Store the weight of each (component, observation) combination
    this->responsibility_matrix = new double*[this->n_obs];
    for (int i = 0; i < this->n_obs; ++i){
        this->responsibility_matrix[i] = new double[this->n_components];
    }

}

void mixtureModel::free_responsibility_matrix(){
    if (this->n_obs != -1){
        for (int i = 0; i < this->n_obs; ++i){
            delete this->responsibility_matrix[i];
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
}

mixtureModel::mixtureModel(const mixtureModel& m){
    this->n_components = m.n_components;
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
    
    this->shared_params_groups = m.shared_params_groups;
    this->shared_params_callbacks = m.shared_params_callbacks;

    for (int i = 0; i < m.n_components; ++i){
        this->weights.push_back(m.weights[i]);
        this->dists.push_back(m.dists[i]);
    }
    this->print_lls = false;
    this->print_dists = false;
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

bool mixtureModel::set_shared_params(vector<int> dist_inds, shared_params_callback callback_func){
    vector<double> meta_params_empty;
    return this->set_shared_params(dist_inds, callback_func, meta_params_empty);
}

bool mixtureModel::set_shared_params(vector<int> dist_inds, shared_params_callback callback_func,
    vector<double> meta_params){
    
    // Ensure that all  dists in the set are valid indices
    set<int> dist_inds_set;
    for (vector<int>::iterator i = dist_inds.begin(); i != dist_inds.end(); ++i){
        if (*i < 0 || *i > this->n_components-1){
            fprintf(stderr, "ERROR: invalid distribution index given in shared parameter set\n");
            return false;
        }
        dist_inds_set.insert(*i);
    }
    this->shared_params_groups.push_back(dist_inds_set);
    this->shared_params_callbacks.push_back(callback_func);
    this->shared_params.push_back(meta_params);
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

    // How many observations in data set?
    long int n_sites = obs.size();
    
    if (n_sites == 0){
        fprintf(stderr, "ERROR: cannot fit mixture model on 0 sites\n");
        exit(1);
    }

    // How many component distributions?
    int n_components = this->n_components;
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
    double** mean_sums = new double*[n_components];
    double member_weight_sums[n_components];
    for (int j = 0; j < n_components; ++j){
        mean_sums[j] = new double[obs[0].size()];
    }
    
    // How small can the variance be before the component is zeroed out? 
    float var_thresh = 0.0;

    // Has the model reached a point where it only contains one component?
    bool one_component = false;
    
    // Make observation weights sum to 1
    double obs_weight_sum = 0.0; 
    for (int i = 0; i < obs_weights.size(); ++i){
        obs_weight_sum += obs_weights[i];
    }
    // Try to prevent underflow by scaling things up by this number
    // (will divide later)
    double obs_weight_scale = (double)obs_weights.size();
    for (int i = 0; i < obs_weights.size(); ++i){
        //obs_weights[i] /= obs_weight_sum;
        obs_weights[i] = pow(2, log2(obs_weights[i]) - log2(obs_weight_sum) + log2(obs_weight_scale));
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
                //member_weight_sums[j] += this->responsibility_matrix[i][j];
                for (int k = 0; k < obs[i].size(); ++k){
                    mean_sums[j][k] += pow(2, log2(this->responsibility_matrix[i][j]) + 
                        log2(obs[i][k]) + log2(obs_weights[i]));
                    //mean_sums[j][k] += this->responsibility_matrix[i][j] * obs[i][k];
                    //mean_sums[j][k] += this->responsibility_matrix[i][j] * obs[i][k] * obs_weights[i];
                }
            }
        }
        
        /* 
        for (int i = 0; i < n_sites; ++i){
            for (int j = 0; j < n_components; ++j){
                if (j > 0){
                    fprintf(stdout, "\t");
                }
                fprintf(stdout, "%f", this->responsibility_matrix[i][j]);
            }
            fprintf(stdout, "\n");
        }
        exit(0);
        */

        // M-step: update parameters
        
        for (int j = 0; j < n_components; ++j){
            if (weights[j] > 0){
                // Distribution weight
                //double new_weight = member_weight_sums[j] / (double)n_sites;
                // This is already normalized by including the (normalized) observation weights
                // in the sum
                double new_weight = member_weight_sums[j] / obs_weight_scale;
                if (isnan(new_weight)){
                    new_weight = 0.0;
                }
                this->weights[j] = new_weight;
            
                if (new_weight > 0){
                    bool can_update = true;

                    // Compute summary statistics (mean and variance of each dimension of observations)
                    vector<double> means;
                    vector<double> vars;
                    for (int k = 0; k < obs[0].size(); ++k){
                        double mean = mean_sums[j][k] / member_weight_sums[j];
                        double var = 0.0;
                        for (int i = 0; i < obs.size(); ++i){
                            //var += this->responsibility_matrix[i][j] * pow(obs[i][k] - mean, 2);
                            //var += pow(log2(this->responsibility_matrix[i][j]) + 
                            //    2*log2(obs[i][k] - mean) + 
                            //    log2(obs_weights[i]), 2);
                            
                            double varlog = log2(this->responsibility_matrix[i][j]) + 
                                2*log2(abs(obs[i][k] - mean)) + log2(obs_weights[i]);
                            var += pow(2, varlog);
                            
                            //var += this->responsibility_matrix[i][j] * pow(obs[i][k] - mean, 2) * obs_weights[i];
                            
                            //var += pow(2, log2(this->responsibility_matrix[i][j]) + 
                            //    2*log2(obs[i][k] - mean) + 
                            //    log2(obs_weights[i]));
                        }
                        var /= weights[j];
                        var /= obs_weight_scale;
                        means.push_back(mean);
                        vars.push_back(var);
                        if (var <= 1e-8){
                            can_update = false;
                        }
                    }
                    
                    if (can_update){
                        // Allow distribution to update using method of moments
                        can_update = this->dists[j].update( means, vars );
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
        
        // Now that we have updated weights and parameters,
        // handle all shared parameters - let whatever external
        // function is set reconcile the parameters of the distrinbutions
        // in the group. 
        for (int group = 0; group < this->shared_params_groups.size(); ++group){

            vector<mixtureDist*> dists_in_group;
            vector<double> dists_in_group_weights;
            double weightsum = 0.0;
            for (set<int>::iterator i = this->shared_params_groups[group].begin(); 
                i != this->shared_params_groups[group].end(); ++i){    
                
                dists_in_group.push_back(&this->dists[*i]);
                dists_in_group_weights.push_back(this->weights[*i]);
                weightsum += this->weights[*i];
            }
            for (int i = 0; i < dists_in_group_weights.size(); ++i){
                dists_in_group_weights[i] /= weightsum;
            }
            this->shared_params_callbacks[group](dists_in_group, 
                dists_in_group_weights, this->shared_params[group]);
        }
        
        for (int j = 0; j < this->weights.size(); ++j){
            if (this->weights[j] == 1.0){
                one_component = true;
                break;
            }
        }
        
        double loglik = this->compute_loglik(obs); 

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
        delta = loglik-loglik_prev;
        if (isnan(delta)){
            fprintf(stderr, "delta NAN %f %f\n", loglik_prev, loglik);
            this->print();
            exit(1);
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
        delete mean_sums[j];
    }
    delete[] mean_sums;
    
    // Store log likelihood of fit model
    this->loglik = loglik_prev;
    
    // Compute BIC and AIC
    this->compute_bic(obs.size() * obs[0].size());    
    return loglik_prev;
}

double mixtureModel::compute_loglik(const vector<vector<double> >& obs){
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
            row_sum += pow(2, ll_row[j] - ll_row_max);
        }
        loglik += log2(row_sum) + ll_row_max;
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
    
    double obs_weight_sum = 0.0;
    for (int i = 0; i < obs_weights.size(); ++i){
        obs_weight_sum += obs_weights[i];
    }
    double obs_weight_scale = (double)obs_weights.size();
    for (int i = 0; i < obs_weights.size(); ++i){
        obs_weights[i] = pow(2, log2(obs_weights[i]) - log2(obs_weight_sum) + log2(obs_weight_scale));
    }

    for (int j = 0; j < obs[0].size(); ++j){
        means.push_back(0.0);
        vars.push_back(0.0);
    }
    double denom = (double)obs.size();
    for (int i = 0; i < obs.size(); ++i){
        for (int j = 0; j < obs[i].size(); ++j){
            means[j] += obs_weights[i] * obs[i][j];
        }
    }
    
    for (int j = 0; j < means.size(); ++j){
        means[j] /= obs_weight_scale;
    }
    
    for (int i = 0; i < obs.size(); ++i){
        for (int j = 0; j < obs[i].size(); ++j){
            vars[j] += pow(obs[i][j] - means[j], 2) * obs_weights[i];
        }
    }
    for (int j = 0; j < vars.size(); ++j){
        vars[j] /= obs_weight_scale;
    }

    /*
    // Transform into columns & use Welford's algorithm to get mean & variance in one pass
    vector<vector<float> > cols;
    for (int k = 0; k < obs[0].size(); ++k){
        vector<float> col;
        cols.push_back(col);
    }

    for (int i = 0; i < obs.size(); ++i){
        for (int k = 0; k < obs[i].size(); ++k){
            cols[k].push_back(obs[i][k]);    
        }
    }
    
    vector<double> means;
    vector<double> vars;
    for (int k = 0; k < cols.size(); ++k){
        pair<double, double> mu_var = welford(cols[k]);
        means.push_back(mu_var.first);
        vars.push_back(mu_var.second);
    }
    */

    // Fit the model (only one distribution)
    this->dists[0].update( means, vars );
    
    this->weights.clear();
    this->weights.push_back(1.0);

    // Compute the log likelihood of the model
    double loglik = compute_loglik(obs);
    this->loglik = loglik;
    this->compute_bic(obs.size() * obs[0].size());
    return loglik;
}


void mixtureModel::print(){
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
    fprintf(stderr, "log2 likelihood = %.3f BIC %.3f AIC %.3f\n", this->loglik,
        this->bic, this->aic);
    for (int i = 0; i < this->n_components; ++i){
        fprintf(stderr, "   Component weight: %.3f\n", this->weights[i]);
        this->dists[i].print(2);
    }
}

