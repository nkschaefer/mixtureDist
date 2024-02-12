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
std::pair<double, double> welford(std::vector<double>& vals);
std::pair<double, double> welford_weights(std::vector<double>& vals,
    std::vector<double>& weights, bool freq_weights);

// ===== PROBABILITY DISTRIBUTIONS =====

double dpois(double x, double l);
double ppois(double x, double l);
double dgamma(double x, double a, double b);
double pgamma(double x, double a, double b);
std::pair<double, double> gamma_moments(double mu, double var);
double dexp(float x, float l);
double dbetabin(double x, double n, double a, double b);
double dbeta(double x, double a, double b);
double pbeta(double x, double a, double b);
std::pair<double, double> beta_moments(double mean, double var);
double dnbinom(int x, int mu, double phi);
double pnbinom(int x, int mu, double phi);
std::pair<int, double> nbinom_moments(double mean, double var);
double dbinom(double n, double k, double p);
double pbinom(int n, int k, double p);
double dmultinom(const std::vector<double>& x, const std::vector<double>& p);
double ddirichlet(const std::vector<double>& x, const std::vector<double>& alpha);
double dhyper(int N, int K, int n, int k);
double dnorm(double x, double mu, double sigma);
double pnorm(double x, double mu, double sigma);
double dbinom_approx(int n, int k, double frac);
double pchisq(double x, double df);

#endif
