#ifndef _MIXTUREDIST_FUNCTIONS_H
#define _MIXTUREDIST_FUNCTIONS_H
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

double binom_coef_log(double n, double k);
double factorial(int n);
double lbetaf(double a, double b);
long double betaf(double a, double b);
double digamma(double x);
std::pair<double, double> welford(std::vector<double>& vals);
std::pair<double, double> welford_weights(std::vector<double>& vals,
    std::vector<double>& weights, bool freq_weights);
double percentile(std::vector<double>& vec, double perc);
std::pair<double, double> int2gauss(double mu1, double mu2, double s1, double s2);
double ov2gauss(double mu1, double mu2, double s1, double s2);
bool bimod_test(std::vector<double>& vals);

// ===== RANDOM SAMPLE FUNCTIONS =====
double rnbinom(double mu, double phi);
double rexp(double lambda);
void rdirichlet(std::vector<double>& params, std::vector<double>& results);

// ===== PROBABILITY DISTRIBUTIONS =====

double dpois(double x, double l);
double ppois(double x, double l);
double dgamma(double x, double a, double b);
double pgamma(double x, double a, double b);
std::pair<double, double> gamma_moments(double mu, double var);
double dexp(float x, float l);
double pexp(double x, double l);
double dbetabin(double x, double n, double a, double b);
double dbeta(double x, double a, double b);
double pbeta(double x, double a, double b);
std::pair<double, double> beta_moments(double mean, double var);
double dnbinom(double x, double mu, double phi);
double pnbinom(double x, double mu, double phi);
std::pair<double, double> nbinom_moments(double mean, double var);
double dbinom(double n, double k, double p);
double pbinom(int n, int k, double p);
double dmultinom(const std::vector<double>& x, const std::vector<double>& p);
double ddirichlet(const std::vector<double>& x, const std::vector<double>& alpha);
double ddirichletmultinom(const std::vector<double>& x, const std::vector<double>& alpha);
double dhyper(int N, int K, int n, int k);
double dnorm(double x, double mu, double sigma);
double pnorm(double x, double mu, double sigma);
double dbinom_approx(int n, int k, double frac);
double pchisq(double x, double df);

#endif
