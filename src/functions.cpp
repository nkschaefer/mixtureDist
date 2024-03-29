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
    return ((lgammaf(a) + lgammaf(b)) - lgammaf(a + b))/M_LN2;
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
    double term2 = lgammaf(a) / log(2);
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
    int mu = (int)round(mean);
    double phi = pow((double)mean, 2) / (var - mean);
    return make_pair(mu, phi);
}

/**
 * Log2 binomial PDF of k successes in n trials with success probability p
 */
double dbinom(double n, double k, double p){
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
 * Log2 Dirichlet PDF
 */
double ddirichlet(const vector<double>& x, const vector<double>& alpha){
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
        //exit(1);
        return log(0);
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
    return (-0.5*pow(a,2)) / log(2) - log2(sigma) - 
        log2(sqrt(2*3.14159265358979));
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
