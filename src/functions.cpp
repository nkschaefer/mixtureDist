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
#include <float.h>
#include "functions.h"
#include "incbeta/incbeta.h"
#include "cdflib/cdflib.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * NOTE: as a standard convention throughout this file, 
 * all PDFs/PMFs are returned as log2(result)
 * while all CDFs are returned as result.
 */

// ===== General functions (not probability distributions) =====

/**
 * Log2 binomial coefficient (n choose k)
 * Uses Stirling's approximation
 */
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
    

    // https://math.stackexchange.com/questions/64716/approximating-the-logarithm-of-the-binomial-coefficient
    // Use Stirling's approximation
    return nf*log2(nf) - kf*log2(kf) - 
        (nf-kf)*log2(nf-kf) + 0.5*(log2(nf) - log2(kf) - 
        log2(nf - kf) - log2(2*M_PI));
}

/**
 * Log2 factorial of n
 * Uses Stirling's approximation
 */
double factorial(int n){
    if (n == 0 || n == 1){
        return 0;      
    }
    else{
       return 0.5*(log2(2) + log2(M_PI) + log2(n)) + 
           (float)n*(log2(n) - log2(exp(1)));
    }
}

/**
 * Log2 of the beta function
 * Valid for all real numbers
 */
double lbetaf(double a, double b){
    int intptr;
    return ((lgammaf_r(a, &intptr) + lgammaf_r(b, &intptr)) - lgammaf_r(a + b, &intptr))/M_LN2;
}

/**
 * Regular (non-log) beta function
 */
long double betaf(double a, double b){
    long double term1 = tgammal(a) / tgammal(a+b);
    long double term2 = tgammal(b) / tgammal(a+b);
    return term1*term2;
}

/**
 * Digamma function (NOTE: base e, not base 2)
 */
double digamma(double x){
    double xcpy = x;
    return psi(&x);
}

/**
 * Welford's algorithm for computing mean and variance in a single pass
 * Returns pair of (mean, variance)
 */
pair<double, double> welford(vector<double>& vals){
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

/**
 * Welford's algorithm for computing mean and variance in a single pass
 * but with weights included
 *
 * From West (1979) as shown on https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Weighted_incremental_algorithm
 *
 * Returns pair of (mean, variance)
 *
 * Third argument = true if weights represent frequency; false if they represent reliability.
 */
pair<double, double> welford_weights(vector<double>& vals, 
    vector<double>& weights, bool freq_weights){
    double w_sum = 0.0;
    double w_sum2 = 0.0;
    double mean = 0.0;
    double S = 0.0;

    for (int i = 0; i < vals.size(); ++i){
        w_sum = w_sum + weights[i];
        w_sum2 = w_sum2 + weights[i] * weights[i];
        double mean_old = mean;
        mean = mean_old + (weights[i] / w_sum) * (vals[i] - mean_old);
        S = S + weights[i] * (vals[i] - mean_old) * (vals[i] - mean);
    }
    
    double pop_var = S / w_sum;
    double var;
    if (freq_weights){
        var = S / (w_sum - 1.0);
    }
    else{
        var = S / (w_sum - w_sum2 / w_sum);
    }
    return make_pair(mean, var);
}

/**
 * Function to return the value from a vector at the given percentile of the
 * empirical distribution.
 */
double percentile(vector<double>& vec, double quant){
    sort(vec.begin(), vec.end());
    double tot = 0.0;
    for (int i = 0; i < vec.size(); ++i){
        tot += vec[i];
    }
    double runtot = 0.0;
    for (int i = 0; i < vec.size(); ++i){
        runtot += vec[i];
        if (runtot >= tot*quant){
            return vec[i];
        }
    }
    return vec[vec.size()-1];
}

/**
 * Finds the point of intersection between two gaussians.
 * Returns two possible solutions
 */
pair<double, double> int2gauss(double mu1, double mu2, double sigma1, double sigma2){
    if (mu2 < mu1){
        // Flip.
        double tmp = mu2;
        double tmp2 = sigma2;
        mu2 = mu1;
        sigma2 = sigma1;
        mu1 = tmp;
        sigma1 = tmp2;
    }
    double sig21 = sigma1*sigma1;
    double sig22 = sigma2*sigma2;
    double mu21 = mu1*mu1;
    double mu22 = mu2*mu2;
    
    double a = (-1.0/sig21 + 1.0/sig22)/2;
    double b = -mu2/sig22 + mu1/sig21;
    double c = (mu22/sig22 - mu21/sig21 + log(sig22/sig21))/2;
    
    // Solve ax^2 + bx + c
    double sol1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
    double sol2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
    return make_pair(sol1, sol2);
}

/**
 * Determine the percent of overlap between two Gaussians
 */
double ov2gauss(double mu1, double mu2, double sigma1, double sigma2){
    
    if (mu2 < mu1){
        // Flip
        double tmp1 = mu1;
        mu1 = mu2;
        mu2 = tmp1;
        double tmp2 = sigma1;
        sigma1 = sigma2;
        sigma2 = tmp2;
    }
    pair<double, double> intpts = int2gauss(mu1, mu2, sigma1, sigma2);
    double ov1 = (1.0 - pnorm(intpts.first, mu1, sigma1) + pnorm(intpts.first, mu2, sigma2));
    double ov2 = (1.0 - pnorm(intpts.second, mu1, sigma1) + pnorm(intpts.second, mu2, sigma2));
    if (ov1 >= 1.0 || ov1 <= 0.0){
        return ov2;
    }
    else{
        return ov1;
    }
}

/**
 * Determine whether a distribution is bimodal, according to Searle's bimodality
 * test
 */
bool bimod_test(vector<double>& data){
    // First, get mean and variance
    pair<double, double> mu_var = welford(data);
    
    double n = (double)data.size();

    double sigma = sqrt(mu_var.second);

    // Now get skewness and kurtosis
    double skew = 0.0;
    double kurt = 0.0;
    
    double frac = 1.0/(n-1.0);

    for (int i = 0; i < data.size(); ++i){
        double diff = (data[i] - mu_var.first)/sigma;
        double diff3 = pow(diff, 3);
        double diff4 = diff3*diff;
        skew += frac*diff3;
        kurt += frac*diff4;
    }

    // Transform kurtosis to excess kurtosis
    kurt -= 3.0;
    
    double scale = (3*(n-1.0)*(n-1.0))/((n-2.0)*(n-3.0));

    double bimod_coef = (skew*skew + 1.0)/(kurt + scale);
    if (bimod_coef > 5.0/9.0){
        return true;
    }
    else{
        return false;
    }
}

// ===== Random value generating functions =====

/**
 * Generate a random sample from a negative binomial distribution
 */
double rnbinom(double mu, double phi){
    
    if (mu <= 0 || phi <= 0){
        return -1;
    }
    while(true){
        double rn = (double)rand() / RAND_MAX;
        
        // Use inverse CDF of Negative Binomial
        
        // First convert mu/phi parameterization to r/p
        double r = phi;
        double p = phi/(mu+phi);
        
        // Treat as inv CDF
        int which = 2;
        // CDF value
        double P = rn;
        // 1 - CDF value
        double Q = 1.0-rn;
        // What they call S wikipedia calls r
        // What they call PR wikipedia calls p
        // This will contain 1-p.
        double OMPR = 1-p;
        double result;
        int status;
        double bound;
        
        cdfnbn(&which, &P, &Q, &result, &r, &p, &OMPR, &status, &bound); 
        if (status != 0){
            //return -1; 
        }
        else{
            return result;    
        }
    }
}

/**
 * Generate a random sample from an exponential distribution
 */
double rexp(double lambda){
    double r = (double)rand() / RAND_MAX;
    return (-log(1.0 - r)/lambda);

}

// ===== Probability distributions and related functions =====

/**
 *  Log2 PMF of Poisson distribution at x with mean l
 */
double dpois(double x, double l){
    if (l < 1){
        return -1;
    }
    double xfac = factorial(x);
    return x * log2(l) + -((float)l/log(2)) - xfac;
}

/**
 * CDF of Poisson distribution at x with mean l
 */
double ppois(double x, double l){
    double p;
    double q;
    double s = x;
    double xlam = l;
    cumpoi(&s, &xlam, &p, &q);
    return p;
}

/**
 * Log2 PDF of gamma distribution at x with parameters a (shape) and b (rate)
 */
double dgamma(double x, double a, double b){
    double term1 = a * log2(b);
    int intptr;
    double term2 = lgammaf_r(a, &intptr) / log(2);
    double term3 = (a - 1.0) * log2(x);
    double term4 = (-b * x) * log2(exp(1));
    return (term1-term2) + term3 + term4;
}

/**
 * CDF of gamma distribution
 */
double pgamma(double x, double a, double b){
    // Convert from shape/rate to shape/scale parameterization
    double shape = a;
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
 * Fit gamma distribution (shape and rate parameters) using
 * method of moments, given mean and variance
 */
pair<double, double> gamma_moments(double mu, double var){
    double alpha = pow(mu, 2) / var;
    double beta = mu / var;
    return make_pair(alpha, beta);
}

/**
 * Log2 PDF of exponential distribution at x with mean l
 */
double dexp(float x, float l){
    if (l <= 0 || x <= 0){
        return -1;
    }
    return log2(l) - l*x*log2(exp(1));
}

/**
 * CDF of exponential distribution at x with mean l
 */
double pexp(double x, double l){
    return 1.0 - exp(-l*x);
}

/**
 * Log2 PDF of beta-binomial distribution
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
 * Log2 PDF of beta distribution
 */
double dbeta(double x, double a, double b){
    float term1 = (a-1) * log2(x);
    float term2 = (b-1) * log2(1.0-x);
    float term3 = lbetaf(a, b);
    return (term1 + term2) - term3;
}

/** 
 * CDF of beta distribution
 */
double pbeta(double x, double a, double b){
    return incbeta(a, b, x); 
}

/**
 * Fit a beta distribution by method of moments
 */
pair<double, double> beta_moments(double mean, double var){
    double alpha = mean*((mean*(1-mean))/var - 1);
    double beta = (alpha * (1-mean))/mean;
    return make_pair(alpha, beta);
}

/**
 * Log2 negative binomial PDF of r successes in k trials with success probability p
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
 * Wikipedia definition:
 * p = mean/variance
 * r = (mean*p)/(1-p)
 *
 */
double dnbinom(int x, int mu, double phi){
    double term1 = binom_coef_log(x + phi - 1, x);
    double term2 = x * (log2(mu) - log2(mu + phi));
    double term3 = phi * (log2(phi) - log2(mu + phi));
    return term1 + term2 + term3;
}

/**
 * Negative binomial CDF 
 * uses mu/phi parameterization
 */
double pnbinom(int x, int mu, double phi){
    double p = phi/((double)mu+phi);
    // incbeta syntax is a, b, x for I_x(a, b)
    return incbeta(phi, x+1, p); 
}

/**
 * Fit negative binomial distribution using method of moments
 */
pair<int, double> nbinom_moments(double mean, double var){
    if (var < mean){
        // Problem
        var = mean + 1.0;
    }
    int mu = (int)round(mean);
    double phi = pow((double)mean, 2) / (var - mean);
    return make_pair(mu, phi);
}

/**
 * Log2 binomial PDF of k successes in n trials with success probability p
 */
double dbinom(double n, double k, double p){
    if (p <= 0){
        p = DBL_MIN*1e6;
    }
    else if (p >= 1){
        p = 1.0-DBL_MIN*1e6;
    }
    double nchoosek = binom_coef_log(n, k);
    double res =  nchoosek + k*log2(p) + (n-k)*log2(1-p);
    return res;
}

/**
 * CDF of Binomial distribution
 */
double pbinom(int n, int k, double p){
    // Regularized incomplete beta function
    // = B(x,a,b)
    // x = q = 1 - p
    // a = n-k
    // b = 1 + k
    return incbeta(n-k, 1+k,1-p);
}

/**
 * Log2 multinomial PDF
 *
 * Takes a vector of double as input, and a vector of double
 * of corresponding length as parameters
 */
double dmultinom(const vector<double>& x, const vector<double>& p){
    if (x.size() != p.size()){
        fprintf(stderr, "ERROR: dmultinom: length of x != length of p\n");
        exit(1);
    }
    double xsum = 1;
    double term2 = 0;
    double term3 = 0;
    double psum = 0.0;
    for (int i = 0; i < x.size(); ++i){
        double this_p = p[i];
        if (this_p <= 0){
            this_p = DBL_MIN*1e6;
        }
        else if (this_p >= 1){
            this_p = 1.0 - DBL_MIN*1e6;
        }
        /*
        if (p[i] <= 0 || p[i] >= 1){
            fprintf(stderr, "ERROR: dmultinom: p[%d] out of range\n", i);
            for (int j = 0; j < x.size(); ++j){
                fprintf(stderr, "p[%d] = %f\n", j, p[j]);
            }
            exit(1);
        }
        */
        psum += this_p;
        xsum += x[i];
        int intptr;
        term2 += lgammaf_r(x[i] + 1, &intptr);
        term3 += x[i] * log(this_p);
    }
    double term1 = lgammaf(xsum);
    if (isinf(term1) || isnan(term1) || isinf(term2) || isnan(term2) || isinf(term3) || isnan(term3)){
        fprintf(stderr, "term1 %f term2 %f term3 %f\n", term1, term2, term3);
        for (int i = 0; i < x.size(); ++i){
            fprintf(stderr, "%d) x = %f p = %f\n", i, x[i], p[i]);
        }
        exit(1);
    }
    return (term1 - term2 + term3)/log(2);
}

/**
 * Log2 Dirichlet PDF
 */
double ddirichlet(const vector<double>& x, const vector<double>& alpha){
    if (x.size() != alpha.size()){
        fprintf(stderr, "ERROR: ddirichlet: length of x != length of alpha\n");
        exit(1);
    }
    double betaprod = 0.0;
    int intptr;
    double term2 = 0.0;
    double alpha_0 = 0.0;
    double xsum = 0.0;
    for (int i = 0; i < alpha.size(); ++i){
        betaprod += lgammaf_r(alpha[i], &intptr);
        alpha_0 += alpha[i];
        xsum += x[i];
        term2 += (alpha[i]-1) * log(x[i]);
    }
    //if (xsum != 1.0){
    //    fprintf(stderr, "ERROR: ddirichlet: x vector does not sum to 1 (%f)\n", xsum);
    //    exit(1);
    //}
    double betadenom = lgammaf_r(alpha_0, &intptr);
    double beta = betaprod - betadenom;
    return (term2 - beta)/log(2);
}

/**
 * Log2 Dirichlet-multinomial PDF
 *
 * Takes a vector of double as input, and a vector of double
 * of corresponding length as parameters
 */
double ddirichletmultinom(const vector<double>& x, const vector<double>& alpha){
    if (x.size() != alpha.size()){
        fprintf(stderr, "ERROR: ddirichletmultinom: length of x != length of alpha\n");
        exit(1);
    }
    
    double sum = 0.0;
    double xtot = 0.0;
    double a_0 = 0.0;
    int intptr;
    for (int i = 0; i < x.size(); ++i){
        double term1 = lgammaf_r(x[i] + alpha[i], &intptr);
        double term2 = lgammaf_r(alpha[i], &intptr);
        double term3 = lgammaf_r(x[i] + 1, &intptr);
        sum += term1 - term2 - term3;
        xtot += x[i];
        a_0 += alpha[i];
    } 

    double term1a = lgammaf_r(a_0, &intptr);
    double term2a = lgammaf_r(xtot+1, &intptr);
    double term3a = lgammaf_r(xtot+a_0, &intptr);
    sum += term1a + term2a - term3a;
    return sum/log(2);
}

/** 
 * Log2 PMF of hypergeometric distribution
 * N population size
 * K success states in population
 * n sample size
 * k sample successes
 */
double dhyper(int N, int K, int n, int k){
    return (binom_coef_log(K, k) + 
        binom_coef_log(N-K, n-k)) - binom_coef_log(N,n);    
}

/**
 * Log2 PDF of normal distribution
 */
double dnorm(double x, double mu, double sigma){
    double a = (x - mu) / sigma;
    double ret = (-0.5*pow(a,2)) / log(2) - log2(sigma) - 
        log2(sqrt(2*3.14159265358979));
    if (isinf(ret) || isnan(ret)){
        fprintf(stderr, "ERROR: dnorm: %f %f %f\n", x, mu, sigma);
        exit(1);
    }
    return ret;
}

/**
 * CDF of normal distribution
 */
double pnorm(double x, double mu, double sigma){
    return 0.5 * erfc(-((x-mu)/sigma) * M_SQRT1_2);
}

/**
 * Log2 PDF of Normal approximation to the binomial distribution
 */
double dbinom_approx(int n, int k, double frac){
    double mu = (double)n * frac;
    double sigma = sqrt((double)n * frac * (1.0-frac));
    return dnorm(k, mu, sigma);
}

/**
 * CDF of Chi-Squared distribution with df degrees of freedom
 */
double pchisq(double x, double df){
    double p;
    double q;
    cumchi(&x, &df, &p, &q);
    return p;
}
