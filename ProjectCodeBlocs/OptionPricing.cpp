
/*
	Paul Wattellier, Mira Maamari, Denisa Draghian, Kacim Younsi
	27/01/2022
	Assignment C++ - End of term Project
	OptionPricing.cpp
*/


#include <cmath>
#include <math.h>
#include <iostream>
#include <random>
std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0,1.0);
#include "OptionPricing.h"



double norm_cdf(double x) { return 0.5 * (1.0 + erf(x / sqrt(2.0))); }

double norm_pdf(double x) {
    return (1.0 / sqrt(2 * M_PI)) * exp(-0.5 * x * x);
}

double rand_normal() {
    return distribution(generator);}


//*********************
//*********************
//	Option Parent Class
//*********************
//*********************

// Constructor
Option::Option(double S, double K, double r, double sigma, double T) :
    S(S), K(K), r(r), sigma(sigma), T(T) {}

Option::~Option() {}
// Getter methods
double Option::getS() { return S; }
double Option::getK() { return K; }
double Option::getR() { return r; }
double Option::getSigma() { return sigma; }
double Option::getT() { return T; }

// Setter methods
void Option::setS(double S) { this->S = S; }
void Option::setK(double K) { this->K = K; }
void Option::setR(double r) { this->r = r; }
void Option::setSigma(double sigma) { this->sigma = sigma; }
void Option::setT(double T) { this->T = T; }

// Calculate d1 and d2
double Option::d1() { return (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T)); }
double Option::d2() { return d1() - sigma * sqrt(T); }

// Function to calculate the delta of the option
double Option::delta() { return norm_cdf(d1()); }

// Function to calculate the gamma of the option
double Option::gamma() { return norm_pdf(d1()) / (S * sigma * sqrt(T)); }

// Function to calculate the theta of the option
double Option::theta() {
    return (-S * norm_pdf(d1()) * sigma) / (2 * sqrt(T)) - (r * K * exp(-r * T) * norm_cdf(d2()));
}

// Function to calculate the vega of the option
double Option::vega() { return S * sqrt(T) * norm_pdf(d1()); }

// Function to calculate the rho of the option
double Option::rho() { return K * T * exp(-r * T) * norm_cdf(d2()); }


//**********************************
//	European Option
//**********************************

EuropeanOption::EuropeanOption(double S, double K, double r, double sigma, double T) :
    Option(S, K, r, sigma, T) {}
EuropeanOption::~EuropeanOption() {}
//**********************************
//	European Call
//**********************************



// Constructor
EuropeanCall::EuropeanCall(double S, double K, double r, double sigma, double T) :
    EuropeanOption(S, K, r, sigma, T) {}
EuropeanCall::~EuropeanCall() {}
// Price a European call option
double EuropeanCall::price() { return getS() * norm_cdf(d1()) - getK() * exp(-getR() * getT()) * norm_cdf(d2()); }

// Replication strategy
void EuropeanCall::replicate() {
    double num_shares = norm_cdf(d1()); // Number of shares of the underlying asset
    double num_bonds = exp(-getR() * getT()) * norm_cdf(d2()); // Number of risk-free bonds
    std::cout << "To replicate a European call option, hold " << num_shares << " shares of the underlying asset and " << num_bonds << " risk-free bonds." << std::endl;
    }


//**********************************
//	European Put
//**********************************

// Constructor
EuropeanPut::EuropeanPut(double S, double K, double r, double sigma, double T) :
        EuropeanOption(S, K, r, sigma, T) {}
EuropeanPut::~EuropeanPut() {}
// Price a European put option
double EuropeanPut::price() { return getK() * exp(-getR() * getT()) * norm_cdf(-d2()) - getS() * norm_cdf(-d1()); }

void EuropeanPut::replicate() {
    double num_shares = -norm_cdf(-d1()); // Number of shares of the underlying asset
    double num_bonds = exp(-getR() * getT()) * norm_cdf(-d2()); // Number of risk-free bonds
    std::cout << "To replicate a European put option, hold " << num_shares << " shares of the underlying asset and " << num_bonds << " risk-free bonds." << std::endl;}

//*******************************************
//*******************************************
//	European Option with Stochastic Volatility
//*******************************************
//*******************************************


StochasticEuropeanOption::StochasticEuropeanOption(double S, double K, double r, double v0, double k, double theta,double sigma_v, double T) :
        Option(S, K, r, 0, T), v0(v0), k(k), theta(theta), sigma_v(sigma_v), n(n) {dt = T/n;}
StochasticEuropeanOption::~StochasticEuropeanOption() {}
// Getter methods
double StochasticEuropeanOption::getV0() { return v0; }
double StochasticEuropeanOption::getK() { return k; }
double StochasticEuropeanOption::getTheta() { return theta; }
double StochasticEuropeanOption::getSigma() { return sigma_v; }
double StochasticEuropeanOption::getn() { return n; }
double StochasticEuropeanOption::getdt() { return dt; }

// Setter methods
void StochasticEuropeanOption::setV0(double v0) { this->v0 = v0; }
void StochasticEuropeanOption::setK(double k) { this->k = k; }
void StochasticEuropeanOption::setTheta(double theta) { this->theta = theta; }
void StochasticEuropeanOption::setSigma(double sigma_v) { this->sigma_v = sigma_v; }
void StochasticEuropeanOption::setn(double n) { this->n = n; }


//**********************************
//	European Call with Stochastic Volatility
//**********************************

// Constructor
StochasticEuropeanCall::StochasticEuropeanCall(double S, double K, double r, double v0, double k, double theta, double sigma, double T) :
    StochasticEuropeanOption(S, K, r, v0, k, theta, sigma, T) {}
StochasticEuropeanCall::~StochasticEuropeanCall() {}
// Function to price a European call option with stochastic volatility
double StochasticEuropeanCall::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double v_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            v_t = getV0();
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * v_t) * getdt() + sqrt(v_t) * rand_normal());
                v_t = v_t + getK() * (getTheta() - v_t) * getdt() + getSigma() * rand_normal();
            }
            payoff = std::max(S_t - getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }

//**********************************
//	European Put with Stochastic Volatility
//**********************************

// Constructor
StochasticEuropeanPut::StochasticEuropeanPut(double S, double K, double r, double v0, double k, double theta, double sigma, double T) :
    StochasticEuropeanOption(S, K, r, v0, k, theta, sigma, T) {}

StochasticEuropeanPut::~StochasticEuropeanPut() {}

// Function to price a European put option with stochastic volatility
double StochasticEuropeanPut::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double v_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            v_t = getV0();
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * v_t) * getdt() + sqrt(v_t) * rand_normal());
                v_t = v_t + getK() * (getTheta() - v_t) * getdt() + getSigma() * rand_normal();
            }
            payoff = std::max(getK() - S_t, 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }



//**********************************
//**********************************
//	Asian Option
//**********************************
//**********************************


AsianOption::AsianOption(double S, double K, double r, double sigma, double T, int n) :
        Option(S, K, r, sigma, T), n(n) {
            dt = T / n;}
AsianOption::~AsianOption() {}
double AsianOption::getn() { return n; }
double AsianOption::getdt() { return dt; }


//**********************************
//	Asian Call Option
//**********************************

AsianCall::AsianCall(double S, double K, double r, double sigma, double T, int n) :
    AsianOption(S, K, r, sigma, T, n) {}
AsianCall::~AsianCall() {}
// Function to price an Asian call option using Monte Carlo methods
double AsianCall::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double average;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            average = 0;
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * sqrt(getdt()) * rand_normal());
                average += S_t;
            }
            average /= getn();
            payoff = std::max(average - getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }


//**********************************
//	Asian Put Option
//**********************************

AsianPut::AsianPut(double S, double K, double r, double sigma, double T, int n) :
    AsianOption(S, K, r, sigma, T, n) {}
AsianPut::~AsianPut() {}
// Function to price an Asian put option using Monte Carlo methods
double AsianPut::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double average;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            average = 0;
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * sqrt(getdt()) * rand_normal());
                average += S_t;
            }
            average /= getn();
            payoff = std::max(getK() - average, 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }


//**********************************
//**********************************
//	Look-back Option
//**********************************
//**********************************

// Constructor
LookbackOption::LookbackOption(double S, double K, double r, double sigma, double T, int n) :
    Option(S, K, r, sigma, T), n(n) {dt = T/n; }
LookbackOption::~LookbackOption() {}
// Getter methods
double LookbackOption::getn() {return n;}
double LookbackOption::getdt() {return dt;}

// Setter methods
void LookbackOption::setn(double n) { this->n = n; }
void LookbackOption::setdt(double dt) { this->dt = dt; }

//**********************************
//	Look back Put Option
//**********************************

// Constructor
LookbackPut::LookbackPut(double S, double K, double r, double sigma, double T, int n) :
    LookbackOption(S, K, r, sigma, T, n) {}
LookbackPut::~LookbackPut() {}
// Function to price a lookback put option using Monte Carlo methods
double LookbackPut::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            double Min_S = getS();
            double Max_S = getS();
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * rand_normal());
                Min_S = std::min(Min_S, S_t);
                Max_S = std::max(Max_S, S_t);
            }
            payoff = std::max(Max_S - S_t, 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }


//**********************************
//	Look back Call Option
//**********************************

// Constructor
LookbackCall::LookbackCall(double S, double K, double r, double sigma, double T, int n) :
    LookbackOption(S, K, r, sigma, T, n) {}

// Function to price a lookback call option using Monte Carlo methods
double LookbackCall::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            double Min_S = getS();
            double Max_S = getS();
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * rand_normal());
                Min_S = std::min(Min_S, S_t);
                Max_S = std::max(Max_S, S_t);
            }
            payoff = std::max(S_t - Min_S, 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }

