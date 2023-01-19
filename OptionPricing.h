
#define OptionPricing_h


#include <cmath>
#include <math.h>
#include <iostream>
#include <random>


//*********************
//*********************
//	Option Parent Class
//*********************
//*********************


class Option
{
private:
    double S;      // Current price of the underlying asset
    double K;      // Strike price
    double r;      // Risk-free interest rate
    double sigma;  // Volatility of the underlying asset
    double T;      // Time to expiration (in years)

public:
    // Constructor
    Option(double S, double K, double r, double sigma, double T);
        virtual ~Option();

    // Getter methods
    double getS();
    double getK();
    double getR();
    double getSigma();
    double getT();

    // Setter methods
    void setS(double S);
    void setK(double K);
    void setR(double r);
    void setSigma(double sigma);
    void setT(double T);

    // Calculate d1 and d2
    double d1();
    double d2();

    // Function to calculate the delta of the option
    double delta();
    // Function to calculate the gamma of the option
    double gamma();
    // Function to calculate the theta of the option
    double theta();
    // Function to calculate the vega of the option
    double vega();
    // Function to calculate the rho of the option
    double rho();
};


//**********************************
//	European Call
//**********************************

class EuropeanCall : public Option
{
public:
    // Constructor
    EuropeanCall(double S, double K, double r, double sigma, double T);
        ~EuropeanCall();

    // Price a European call option
    double price();

    // Replication strategy
    void replicate();
    };

//**********************************
//	European Put
//**********************************
class EuropeanPut : public Option
{
public:
    // Constructor
    EuropeanPut(double S, double K, double r, double sigma, double T);
        ~EuropeanPut();

    // Price a European put option
    double price();

    // Replication strategy
    void replicate();
};


//*******************************************
//*******************************************
//	European Option with Stochastic Volatility
//*******************************************
//*******************************************


class StochasticEuropeanOption : public Option
{
private:
    double v0; // Initial volatility
    double k; // Mean reversion speed of volatility
    double theta; // Long-term mean of volatility
    double rho; // Correlation between asset and volatility
    double sigma; // Volatility of volatility
    double n;
    double dt;

public:
    // Constructor
    StochasticEuropeanOption(double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T);
        ~StochasticEuropeanOption();

    // Getter methods
    double getV0();
    double getK();
    double getTheta();
    double getRho();
    double getSigma();
    double getn();
    double getdt();
    // Setter methods
    void setV0(double v0);
    void setK(double k);
    void setTheta(double theta) ;
    void setRho(double rho);
    void setSigma(double sigma);
    void setn(double n);

};
//**********************************
//	European Call with Stochastic Volatility
//**********************************


class StochasticEuropeanCall : public StochasticEuropeanOption
{
public:
    // Constructor
    StochasticEuropeanCall(double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T);
        ~StochasticEuropeanCall();


    // Function to price a European call option with stochastic volatility
    double price(int num_simulations);
};

//**********************************
//	European Put with Stochastic Volatility
//**********************************

class StochasticEuropeanPut : public StochasticEuropeanOption
{
public:
    // Constructor
    StochasticEuropeanPut(double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T);
        ~StochasticEuropeanPut();

    // Function to price a European put option with stochastic volatility
    double price(int num_simulations);
};


//**********************************
//**********************************
//	Asian Option
//**********************************
//**********************************



class AsianOption : public Option
{
private:
    int n;  // Number of time steps
    double dt; // Time step size
    double payoff; // Payoff of the option
public:
    // Constructor
    AsianOption(double S, double K, double r, double sigma, double T, int n);
        ~AsianOption();
    double getn();
    double getdt();
};

//**********************************
//	Asian Call Option
//**********************************

class AsianCall : public AsianOption
{
public:
    // Constructor
    AsianCall(double S, double K, double r, double sigma, double T, int n);
        ~AsianCall();
    // Function to price an Asian call option using Monte Carlo methods
    double price(int num_simulations);
};

//**********************************
//	Asian Put Option
//**********************************

class AsianPut : public AsianOption
{
public:
    // Constructor
    AsianPut(double S, double K, double r, double sigma, double T, int n);
        ~AsianPut();

    // Function to price an Asian put option using Monte Carlo methods
    double price(int num_simulations);
};



//**********************************
//**********************************
//	Look-back Option
//**********************************
//**********************************

class LookbackOption : public Option
{
private:
    double n;
    double dt;

public:
    // Constructor
    LookbackOption(double S, double K, double r, double sigma, double T, int n);
        ~LookbackOption();

    // Getter methods
    double getn();
    double getdt();

    // Setter methods
    void setn(double n);
    void setdt(double dt);
};

//**********************************
//	Look back Put Option
//**********************************

class LookbackPut : public LookbackOption
{
public:
    // Constructor
    LookbackPut(double S, double K, double r, double sigma, double T, int n);
        ~LookbackPut();

    // Function to price a lookback put option using Monte Carlo methods
    double price(int num_simulations);
};

//**********************************
//	Look back Call Option
//**********************************
class LookbackCall : public LookbackOption
{
public:
    // Constructor
    LookbackCall(double S, double K, double r, double sigma, double T, int n);

    // Function to price a lookback call option using Monte Carlo methods
    double price(int num_simulations);

};

