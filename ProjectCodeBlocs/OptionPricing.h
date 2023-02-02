/*
	Paul Wattellier, Mira Maamari, Denisa Draghia, Kacim Younsi
	27/01/2022
	Assignment C++ - End of term Project
	OptionPricing.h
*/

#define OptionPricing_h
#include <cmath>
#include <math.h>
#include <iostream>
#include <random>
#include <vector>
#include <cstdlib>

//*********************
//*********************
//	Option Parent Class

//*********************
//*********************

class Option
{
protected:
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


    // The Greeks

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

    // Price a European call option using BS formula
    double price();

    // Price a European put option using Monte Carlo simulation
    double price_MonteCarlo(int);

    //Price difference between the BS price and the Monte Carlo one for n simulations
    double difference(int);


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

    // Price a European put option using BS formula
    double price();

    // Price a European put option using Monte Carlo simulation
    double price_MonteCarlo(int);

    //Price difference between the BS price and the Monte Carlo one for n simulations
    double difference(int);

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
    double getk();
    double getTheta();
    double getRho();
    double getSigma();
    double getn();
    double getdt();
    // Setter methods
    void setV0(double v0);
    void setk(double k);
    void setTheta(double theta) ;
    void setRho(double rho);
    void setSigma(double sigma);
    void setn(double n);

};
//**********************************
//	European Call with Stochastic Volatility
//**********************************


class StochasticEuropeanCall_volatility : public StochasticEuropeanOption
{
public:
    // Constructor
    StochasticEuropeanCall_volatility (double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T);
        ~StochasticEuropeanCall_volatility ();


    // Function to price a European call option with stochastic volatility
    double price(int num_simulations);
};

//**********************************
//	European Put with Stochastic Interest rate
//**********************************

class StochasticEuropeanPut_interest  : public StochasticEuropeanOption
{
public:
    // Constructor
    StochasticEuropeanPut_interest (double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T);
        ~StochasticEuropeanPut_interest ();

    // Function to price a European put option with stochastic interest rate
    double price(int num_simulations);
};

//**********************************
//	European Call with Stochastic Interest
//**********************************


class StochasticEuropeanCall_interest : public StochasticEuropeanOption
{
public:
    // Constructor
    StochasticEuropeanCall_interest (double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T);
        ~StochasticEuropeanCall_interest ();


    // Function to price a European call option with stochastic interest rate
    double price(int num_simulations);
};

//**********************************
//	European Put with Stochastic Volatility
//**********************************

class StochasticEuropeanPut_volatility  : public StochasticEuropeanOption
{
public:
    // Constructor
    StochasticEuropeanPut_volatility (double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T);
        ~StochasticEuropeanPut_volatility ();

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
protected:
    int n;  // Number of time steps
    double dt; // Time step size
    double payoff; // Payoff of the option
public:
    // Constructor
    AsianOption(double S, double K, double r, double sigma, double T, int n);
        ~AsianOption();
    // Getter methods
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
protected:
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


//**********************************
//**********************************
//	Binary Option
//**********************************
//**********************************



class BinaryOption : public Option
{
protected:
    int n;  // Number of time steps
    double dt; // Time step size
    double payoff; // Payoff of the option
public:
    // Constructor
    BinaryOption(double S, double K, double r, double sigma, double T, int n);
        ~BinaryOption();
    // Getter methods
    double getn();
    double getdt();
};

//**********************************
//	Binary Call Option
//**********************************

class BinaryCall : public BinaryOption
{
public:
    // Constructor
    BinaryCall(double S, double K, double r, double sigma, double T, int n);
    ~BinaryCall();
    // Function to price a Binary call option using Monte Carlo methods
    double price(int num_simulations);
};

//**********************************
//	Binary Put Option
//**********************************

class BinaryPut : public BinaryOption
{
public:
    // Constructor
    BinaryPut(double S, double K, double r, double sigma, double T, int n);
    ~BinaryPut();
    // Function to price an Binary put option using Monte Carlo methods
    double price(int num_simulations);
};

//**********************************
//**********************************
//	Gap Option
//**********************************
//**********************************

// A gap call option is a European call option that pays off S_t-K1 when S_t>K2.
//Its price is the price of a European Call priced using BS formula + (K2-K1)*exp (-r*T)*N(d2)


//**********************************
//	Gap Call Option
//**********************************
class GapCall : public Option
{private:
    double K1;
public:
    // Constructor
    GapCall(double S, double K, double r, double sigma, double T, double K1);
        ~GapCall();
    // Function to price an Gap call option
    double price();
    //Getter methods
    double getK1();
    };

//**********************************
//	Gap Put Option
//**********************************

// The price of a Gas Put Option is K1*exp(-r*T)*N(-d2)-S_T* N(-d1)
class GapPut : public Option
{private:
    double K1;
public:
    // Constructor
    GapPut(double S, double K, double r, double sigma, double T, double K1);
        ~GapPut();
    // Function to price an Gap call option
    double price();
    //Getter methods
    double getK1();
    };

//**********************************
//**********************************
//	Chooser Option
//**********************************
//**********************************

//The option allows to decide, at a predefined time T1, whether to exercice a call or a put

class ChooserOption : public Option
{
private:
    double T1;
public:
    // Constructor
    ChooserOption(double S, double K, double r, double sigma, double T,double T1);
        ~ChooserOption();
    // Getter methods
    double getT1();
    // Function to price an Chooser option
    double price();
};

//**********************************
//**********************************
//	Barrier Option
//**********************************
//**********************************



class BarrierOption : public Option
{
protected:
    int type;// type=1 if the option is knock-in and type=0 if the option is knock-out
    double n;
    double dt;
    double barrier; // the value of the barrier
public:
    // Constructor
    BarrierOption(double S, double K, double r, double sigma, double T,double n,int type,double barrier);
        ~BarrierOption();
    int getType();
    double getbarrier();
    double getn();
    double getdt();

    // Setter methods
    void setn(double n);
    void setdt(double dt);
    void settype(int type);
    void setbarrier(double);
};


//**********************************
//	Barrier Call Option
//**********************************
class BarrierCall: public BarrierOption
{private:
int type_call_1; // type_call_1=1 if the option is UP-and- and
//type_call_1=0 if the option is DOWN-and-

public:
//Constructor
BarrierCall(double S, double K, double r, double sigma, double T, double n, int type, double barrier, int type_call_1);
~BarrierCall();

//Function to price Barrier Call Option using Monte Carlo simulation
double price(int num_simulations);

//Getter methods
int getType1();

};

//**********************************
//	Barrier Put Option
//**********************************
class BarrierPut: public BarrierOption
{private:
int type_put_1; // type_put_1=1 if the option is UP-and- and
//type_put_1=0 if the option is DOWN-and-

public:
//Constructor
BarrierPut(double S, double K, double r, double sigma, double T, double n, int type, double barrier, int type_call_1);
~BarrierPut();

//Function to price Barrier Put Option using Monte Carlo simulation
double price(int num_simulations);

//Getter methods
int getType1();

};
