/*
	Paul Wattellier, Mira Maamari, Denisa Draghia, Kacim Younsi
	27/01/2022
	Assignment C++ - End of term Project
	OptionPricing.cpp
*/


#include <cmath>
#include <math.h>
#include <iostream>
#include <random>
#include <stdio.h>
#include <vector>
#include <assert.h>
std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0,1.0);
#include "OptionPricing.h"

double PI=3.14159265359;
//Box-Muller algorithm to generate random gaussian numbers
double gaussian_box_muller()
{
  double U1=((double)rand()/(double)RAND_MAX);
  double U2=((double)rand()/(double)RAND_MAX);
  double n=sqrt(-2*log(U1))*cos(2*U2*PI);


  return n;};


double norm_cdf(double x) { return 0.5 * (1.0 + erf(x / sqrt(2.0))); }

double norm_pdf(double x) {
    return (1.0 / sqrt(2 * PI)) * exp(-0.5 * x * x);
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
//	European Call
//**********************************



// Constructor
EuropeanCall::EuropeanCall(double S, double K, double r, double sigma, double T) :
    Option(S, K, r, sigma, T) {}
EuropeanCall::~EuropeanCall() {}

// Price a European call option by BS formula
double EuropeanCall::price() { return getS() * norm_cdf(d1()) - getK() * exp(-getR() * getT()) * norm_cdf(d2()); }

//Price a European call option by Monte Carlo Simulation

double EuropeanCall::price_MonteCarlo(int num_simulations) {
        double sum = 0;
        double S_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < 12; j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getT()/12 + getSigma() * sqrt(getT()/12) * gaussian_box_muller());

            }
            payoff = std::max(S_t - getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }

double EuropeanCall::difference(int num_simulations){

    return abs(price()-price_MonteCarlo(num_simulations));};


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
        Option(S, K, r, sigma, T) {}

EuropeanPut::~EuropeanPut() {}
// Price a European put option
double EuropeanPut::price() { return getK() * exp(-getR() * getT()) * norm_cdf(-d2()) - getS() * norm_cdf(-d1()); }

//Price a European call option by Monte Carlo Simulation

double EuropeanPut::price_MonteCarlo(int num_simulations) {
        double sum = 0;
        double S_t;

        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < 12; j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getT()/12 + getSigma() * sqrt(getT()/12) * gaussian_box_muller());

            }
            payoff = std::max(-S_t + getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }

double EuropeanPut::difference(int num_simulations){
    return abs(price()-price_MonteCarlo(num_simulations));};



void EuropeanPut::replicate() {
    double num_shares = -norm_cdf(-d1()); // Number of shares of the underlying asset
    double num_bonds = exp(-getR() * getT()) * norm_cdf(-d2()); // Number of risk-free bonds
    std::cout << "To replicate a European put option, hold " << num_shares << " shares of the underlying asset and " << num_bonds << " risk-free bonds." << std::endl;}

//*******************************************
//*******************************************
//	European Option with Stochastic Volatility
//*******************************************
//*******************************************


StochasticEuropeanOption::StochasticEuropeanOption(double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T) :
        Option(S, K, r, 0, T), v0(v0), k(k), theta(theta), rho(rho), sigma(sigma), n(n) {dt = T/n;}
StochasticEuropeanOption::~StochasticEuropeanOption() {}
// Getter methods
double StochasticEuropeanOption::getV0() { return v0; }
double StochasticEuropeanOption::getK() { return k; }
double StochasticEuropeanOption::getTheta() { return theta; }
double StochasticEuropeanOption::getRho() { return rho; }
double StochasticEuropeanOption::getSigma() { return sigma; }
double StochasticEuropeanOption::getn() { return n; }
double StochasticEuropeanOption::getdt() { return dt; }

// Setter methods
void StochasticEuropeanOption::setV0(double v0) { this->v0 = v0; }
void StochasticEuropeanOption::setK(double k) { this->k = k; }
void StochasticEuropeanOption::setTheta(double theta) { this->theta = theta; }
void StochasticEuropeanOption::setRho(double rho) { this->rho = rho; }
void StochasticEuropeanOption::setSigma(double sigma) { this->sigma = sigma; }
void StochasticEuropeanOption::setn(double n) { this->n = n; }


//**********************************
//	European Call with Stochastic Volatility
//**********************************

// Constructor
StochasticEuropeanCall_volatility ::StochasticEuropeanCall_volatility (double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T) :
    StochasticEuropeanOption(S, K, r, v0, k, theta, rho, sigma, T) {}
StochasticEuropeanCall_volatility ::~StochasticEuropeanCall_volatility () {}
// Function to price a European call option with stochastic volatility
double StochasticEuropeanCall_volatility ::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double v_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            v_t = getV0();
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * v_t) * getdt() + sqrt(v_t) * gaussian_box_muller());
                v_t = v_t + getK() * (getTheta() - v_t) * getdt() + getSigma() * gaussian_box_muller();
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
StochasticEuropeanPut_volatility ::StochasticEuropeanPut_volatility (double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T) :
    StochasticEuropeanOption(S, K, r, v0, k, theta, rho, sigma, T) {}

StochasticEuropeanPut_volatility ::~StochasticEuropeanPut_volatility () {}

// Function to price a European put option with stochastic volatility
double StochasticEuropeanPut_volatility ::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double v_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            v_t = getV0();
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * v_t) * getdt() + sqrt(v_t) * gaussian_box_muller());
                v_t = v_t + getK() * (getTheta() - v_t) * getdt() + getSigma() * gaussian_box_muller();
            }
            payoff = std::max(getK() - S_t, 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }

//**********************************
//	European Call with Stochastic Interest rate
//**********************************

// Constructor
StochasticEuropeanCall_interest::StochasticEuropeanCall_interest (double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T) :
    StochasticEuropeanOption(S, K, r, v0, k, theta, rho, sigma, T) {}
StochasticEuropeanCall_interest ::~StochasticEuropeanCall_interest () {}
// Function to price a European call option with stochastic interest rate
double StochasticEuropeanCall_interest::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double r_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            r_t = getR();
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getV0()) * getdt() + sqrt(getV0()) * gaussian_box_muller());
                r_t = r_t + getK() * (getTheta() - r_t) * getdt() + getSigma() * gaussian_box_muller();
            }
            payoff = std::max(S_t - getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }

//**********************************
//	European Put with Stochastic Interest rate
//**********************************

// Constructor
StochasticEuropeanPut_interest::StochasticEuropeanPut_interest (double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T) :
    StochasticEuropeanOption(S, K, r, v0, k, theta, rho, sigma, T) {}

StochasticEuropeanPut_interest ::~StochasticEuropeanPut_interest () {}

// Function to price a European put option with stochastic interest rate
double StochasticEuropeanPut_interest ::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double r_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();
            r_t = getR();
            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getV0()) * getdt() + sqrt(getV0()) * gaussian_box_muller());
                r_t = r_t + getK() * (getTheta() - r_t) * getdt() + getSigma() * gaussian_box_muller();
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
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * sqrt(getdt()) * gaussian_box_muller());
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
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * sqrt(getdt()) * gaussian_box_muller());
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
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
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
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
                Min_S = std::min(Min_S, S_t);
                Max_S = std::max(Max_S, S_t);
            }
            payoff = std::max(S_t - Min_S, 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }

//**********************************
//**********************************
//	Binary Option
//**********************************
//**********************************


//Function Heaviside that return 1 if its argument is positive and 0 otherwise

double heaviside(const double& x) {
  if (x >= 0) {
      return 1.0;
  } else {
    return 0.0;
  }
}

//Constructor
BinaryOption::BinaryOption(double S, double K, double r, double sigma, double T, int n) :
        Option(S, K, r, sigma, T), n(n) {
            dt = T / n;}
BinaryOption::~BinaryOption() {}

//Getter methods
double BinaryOption::getn() { return n; }
double BinaryOption::getdt() { return dt; }



//**********************************
//	Binary Call Option
//**********************************

// Constructor
BinaryCall::BinaryCall(double S, double K, double r, double sigma, double T, int n) :
    BinaryOption(S, K, r, sigma, T, n) {}
BinaryCall::~BinaryCall() {}

// Function to price an Binary call option using Monte Carlo methods
double BinaryCall::price(int num_simulations) {
    double sum = 0;
    double S_t;
    double payoff;
    for (int i = 0; i < num_simulations; i++) {
        S_t = getS();
        for (int j = 0; j < 12; j++) {
            S_t = S_t * exp((getR()-0.5 *getSigma() * getSigma())* getT()/12+getSigma()*sqrt(getT()/12) * gaussian_box_muller());
            }

            payoff = heaviside(S_t - getK());
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }


//**********************************
//	Binary Put Option
//**********************************

// Constructor
BinaryPut::BinaryPut(double S, double K, double r, double sigma, double T, int n) :
    BinaryOption(S, K, r, sigma, T, n) {}
BinaryPut::~BinaryPut() {}

// Function to price an Binary call option using Monte Carlo methods
double BinaryPut::price(int num_simulations) {
        double sum = 0;
        double S_t;
        double payoff;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < 12; j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getT()/12 + getSigma() * sqrt(getT()/12) * gaussian_box_muller());

            }
            payoff = heaviside(-S_t + getK());
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }


//**********************************
//**********************************
//	Gap Option
//**********************************
//**********************************


//**********************************
//	Gap Call Option
//**********************************

// Constructor
GapCall::GapCall(double S, double K, double r, double sigma, double T, double K1) :
    Option(S, K, r, sigma, T) {}
GapCall::~GapCall() {}

//Getter methods
double GapCall::getK1(){return K1;};

//Function to price a Gap Call option
double GapCall::price()
{return getS() * norm_cdf(d1()) - getK() * exp(-getR() * getT()) * norm_cdf(d2())+(getK()-getK1())*norm_cdf(d2()); };


//**********************************
//	Gap Put Option
//**********************************

// Constructor
GapPut::GapPut(double S, double K, double r, double sigma, double T, double K1) :
    Option(S, K, r, sigma, T) {}
GapPut::~GapPut() {}

//Getter methods
double GapPut::getK1(){return K1;};
// Function to price a Gap Put option
double GapPut::price() { return getK1()*exp(-getR()*getT())*norm_cdf(-d2())-getS()*norm_cdf(-d1());};



//**********************************
//**********************************
//	Chooser Option
//**********************************
//**********************************


// Constructor
ChooserOption::ChooserOption(double S, double K, double r, double sigma, double T, double T1) :
    Option(S, K, r, sigma, T) {}
ChooserOption::~ChooserOption() {}

// Getter methods
double ChooserOption::getT1(){return T1;};

// Function to price an Chooser option
// The price of a chooser option is the price of an European Call Option with maturity of T-T1+
//the price of an European Put option with a strike of K*exp(-r(T-T1)) and a maturity of T-T1

double ChooserOption::price()
{return EuropeanCall(getS(),getK(),getR(),getSigma(),getT()-getT1()).price()
+EuropeanPut(getS(), getK()*exp(-getR()*(getT()-getT1())),getR(),getSigma(),getT()-getT1()).price();};


//**********************************
//**********************************
//	Barrier Option
//**********************************
//**********************************


BarrierOption::BarrierOption(double S, double K, double r, double sigma, double T, double n, int typ, double barrie) :
    Option(S, K, r, sigma, T), n(n), type(typ), barrier(barrie) {dt=T/n;};

BarrierOption::~BarrierOption() {}
// Getter methods
double BarrierOption::getn() {return n;}
double BarrierOption::getdt() {return dt;}
int BarrierOption::getType() {return type;}
double BarrierOption::getbarrier() {return barrier;}
// Setter methods
void BarrierOption::setn(double n) { this->n = n; }
void BarrierOption::setdt(double dt) { this->dt = dt; }
void BarrierOption::setbarrier(double barrier) {this->barrier=barrier;}

//**********************************
//	Barrier Call Option
//**********************************

//Constructor
BarrierCall::BarrierCall(double S, double K, double r, double sigma, double T, double n,int type,double barrier, int type_call_1) :
    BarrierOption(S, K, r, sigma, T, n,type,barrier) {type_call_1=getType1();}
BarrierCall::~BarrierCall(){}

//Getter methods
int BarrierCall::getType1() {return type_call_1;}


// Function to price a barrier call option using Monte Carlo methods
double BarrierCall::price(int num_simulations)
{if ((getType()==1)&&(getType1()==1)) // The call is up-and-in
{if (getbarrier()<getS()) {std::cout<<"Wrong parameters ";
    return -1;
}
else {double sum=0;
        double S_t;
        double payoff;
        int kick=0;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
                if (S_t>getbarrier()){kick+=1;}
            }
            if (kick==0){payoff =0;}
            else payoff=std::max(S_t-getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }}
if ((getType()==1)&&(getType1()==0)) // The call is down-and-in
{if (getbarrier()>getS()) {std::cout<<"Wrong parameters ";
    return -1;
}
else {double sum=0;
        double S_t;
        double payoff;
        int kick=0;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
                if (S_t<getbarrier()){kick+=1;}
            }
            if (kick==0){payoff =0;}
            else payoff=std::max(S_t-getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }}

if ((getType()==0)&&(getType1()==0)) // The call is down-and-out
{if (getbarrier()>getS()) {std::cout<<"Wrong parameters ";
    return -1;
}
else {double sum=0;
        double S_t;
        double payoff;
        int kick=0;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
                if (S_t<getbarrier()){kick+=1;}
            }
            if (kick==0){payoff =0;}
            else payoff=std::max(S_t-getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }}

if ((getType()==0)&&(getType1()==1)) // The call is up-and-out
{if (getbarrier()<getS()) {std::cout<<"Wrong parameters ";
    return -1;
}
else {double sum=0;
        double S_t;
        double payoff;
        int kick=0;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
                if (S_t>getbarrier()){kick+=1;}
            }
            if (kick>0){payoff =0;}
            else payoff=std::max(S_t-getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());

    }}
    return -1;
};



//**********************************
//	Barrier Put Option
//**********************************


//Constructor
BarrierPut::BarrierPut(double S, double K, double r, double sigma, double T, double n,int type,double barrier, int type_call_1) :
    BarrierOption(S, K, r, sigma, T, n,type,barrier) {type_call_1=getType1();}
BarrierPut::~BarrierPut(){}

//Getter methods
int BarrierPut::getType1() {return type_put_1;}


// Function to price a barrier put option using Monte Carlo methods
double BarrierPut::price(int num_simulations)
{
    if ((getType()==1)&&(getType1()==1)) // The put is up-and-in
{if (getbarrier()<getS()) {std::cout<<"Wrong parameters ";
    return -1;
}
else {double sum=0;
        double S_t;
        double payoff;
        int kick=0;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
                if (S_t>getbarrier()){kick+=1;}
            }
            if (kick==0){payoff =0;}
            else payoff=std::max(-S_t+getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }}
if ((getType()==1)&&(getType1()==0)) // The put is down-and-in
{if (getbarrier()>getS()) {std::cout<<"Wrong parameters ";
    return -1;
}
else {double sum=0;
        double S_t;
        double payoff;
        int kick=0;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
                if (S_t<getbarrier()){kick+=1;}
            }
            if (kick==0){payoff =0;}
            else payoff=std::max(-S_t+getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }}

if ((getType()==0)&&(getType1()==0)) // The put is down-and-out
{if (getbarrier()>getS()) {std::cout<<"Wrong parameters ";
    return -1;
}
else {double sum=0;
        double S_t;
        double payoff;
        int kick=0;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
                if (S_t<getbarrier()){kick+=1;}
            }
            if (kick==0){payoff =0;}
            else payoff=std::max(-S_t+getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }}

if ((getType()==0)&&(getType1()==1)) // The put is up-and-out
{if (getbarrier()<getS()) {std::cout<<"Wrong parameters ";
    return -1;
}
else {double sum=0;
        double S_t;
        double payoff;
        int kick=0;
        for (int i = 0; i < num_simulations; i++) {
            S_t = getS();

            for (int j = 0; j < getn(); j++) {
                S_t = S_t * exp((getR() - 0.5 * getSigma() * getSigma()) * getdt() + getSigma() * gaussian_box_muller());
                if (S_t>getbarrier()){kick+=1;}
            }
            if (kick>0){payoff =0;}
            else payoff=std::max(-S_t+getK(), 0.0);
            sum += payoff;
        }
        return (sum / num_simulations) * exp(-getR() * getT());
    }}
    return -1;
};
















