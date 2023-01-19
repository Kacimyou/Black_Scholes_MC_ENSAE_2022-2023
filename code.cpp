#include <cmath>
#include <math.h>
#include <iostream>
#include <random>

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0,1.0);

double norm_cdf(double x) { return 0.5 * (1.0 + erf(x / sqrt(2.0))); }

double norm_pdf(double x) {
    return (1.0 / sqrt(2 * M_PI)) * exp(-0.5 * x * x);
}

double rand_normal() {
    return distribution(generator);
}


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
    Option(double S, double K, double r, double sigma, double T) :
        S(S), K(K), r(r), sigma(sigma), T(T) {}

    // Getter methods
    double getS() { return S; }
    double getK() { return K; }
    double getR() { return r; }
    double getSigma() { return sigma; }
    double getT() { return T; }

    // Setter methods
    void setS(double S) { this->S = S; }
    void setK(double K) { this->K = K; }
    void setR(double r) { this->r = r; }
    void setSigma(double sigma) { this->sigma = sigma; }
    void setT(double T) { this->T = T; }

    // Calculate d1 and d2
    double d1() { return (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T)); }
    double d2() { return d1() - sigma * sqrt(T); }

    // Function to calculate the delta of the option
    double delta() { return norm_cdf(d1()); }

    // Function to calculate the gamma of the option
    double gamma() { return norm_pdf(d1()) / (S * sigma * sqrt(T)); }

    // Function to calculate the theta of the option
    double theta() {
        return (-S * norm_pdf(d1()) * sigma) / (2 * sqrt(T)) - (r * K * exp(-r * T) * norm_cdf(d2()));
    }

    // Function to calculate the vega of the option
    double vega() { return S * sqrt(T) * norm_pdf(d1()); }

    // Function to calculate the rho of the option
    double rho() { return K * T * exp(-r * T) * norm_cdf(d2()); }

};

class EuropeanCall : public Option
{
public:
    // Constructor
    EuropeanCall(double S, double K, double r, double sigma, double T) :
        Option(S, K, r, sigma, T) {}

    // Price a European call option
    double price() { return getS() * norm_cdf(d1()) - getK() * exp(-getR() * getT()) * norm_cdf(d2()); }

    // Replication strategy
    void replicate() {
        double num_shares = norm_cdf(d1()); // Number of shares of the underlying asset
        double num_bonds = exp(-getR() * getT()) * norm_cdf(d2()); // Number of risk-free bonds
        std::cout << "To replicate a European call option, hold " << num_shares << " shares of the underlying asset and " << num_bonds << " risk-free bonds." << std::endl;
}};

class EuropeanPut : public Option
{
public:
    // Constructor
    EuropeanPut(double S, double K, double r, double sigma, double T) :
        Option(S, K, r, sigma, T) {}

    // Price a European put option
    double price() { return getK() * exp(-getR() * getT()) * norm_cdf(-d2()) - getS() * norm_cdf(-d1()); }

    void replicate() {
        double num_shares = -norm_cdf(-d1()); // Number of shares of the underlying asset
        double num_bonds = exp(-getR() * getT()) * norm_cdf(-d2()); // Number of risk-free bonds
        std::cout << "To replicate a European put option, hold " << num_shares << " shares of the underlying asset and " << num_bonds << " risk-free bonds." << std::endl;}
};

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
    StochasticEuropeanOption(double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T) :
        Option(S, K, r, 0, T), v0(v0), k(k), theta(theta), rho(rho), sigma(sigma), n(n) {dt = T/n;}

    // Getter methods
    double getV0() { return v0; }
    double getK() { return k; }
    double getTheta() { return theta; }
    double getRho() { return rho; }
    double getSigma() { return sigma; }
    double getn() { return n; }
    double getdt() { return dt; }

    // Setter methods
    void setV0(double v0) { this->v0 = v0; }
    void setK(double k) { this->k = k; }
    void setTheta(double theta) { this->theta = theta; }
    void setRho(double rho) { this->rho = rho; }
    void setSigma(double sigma) { this->sigma = sigma; }
    void setn(double n) { this->n = n; }

    // Function to price the option using Monte Carlo methods
    virtual double price(int num_simulations) = 0;
};

class StochasticEuropeanCall : public StochasticEuropeanOption
{
public:
    // Constructor
    StochasticEuropeanCall(double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T) :
        StochasticEuropeanOption(S, K, r, v0, k, theta, rho, sigma, T) {}

    // Function to price a European call option with stochastic volatility
    double price(int num_simulations) {
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
};

class StochasticEuropeanPut : public StochasticEuropeanOption
{
public:
    // Constructor
    StochasticEuropeanPut(double S, double K, double r, double v0, double k, double theta, double rho, double sigma, double T) :
        StochasticEuropeanOption(S, K, r, v0, k, theta, rho, sigma, T) {}

    // Function to price a European put option with stochastic volatility
    double price(int num_simulations) {
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
};






class AsianOption : public Option
{
private:
    int n;  // Number of time steps
    double dt; // Time step size
    double payoff; // Payoff of the option
public:
    // Constructor
    AsianOption(double S, double K, double r, double sigma, double T, int n) :
        Option(S, K, r, sigma, T), n(n) {
            dt = T / n;
        }
    double getn() { return n; }
    double getdt() { return dt; }
};

class AsianCall : public AsianOption
{
public:
    // Constructor
    AsianCall(double S, double K, double r, double sigma, double T, int n) :
        AsianOption(S, K, r, sigma, T, n) {}

    // Function to price an Asian call option using Monte Carlo methods
    double price(int num_simulations) {
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
};

class AsianPut : public AsianOption
{
public:
    // Constructor
    AsianPut(double S, double K, double r, double sigma, double T, int n) :
        AsianOption(S, K, r, sigma, T, n) {}

    // Function to price an Asian put option using Monte Carlo methods
    double price(int num_simulations) {
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
};

class LookbackOption : public Option
{
private:
    double n;
    double dt;

public:
    // Constructor
    LookbackOption(double S, double K, double r, double sigma, double T, int n) :
        Option(S, K, r, sigma, T), n(n) {dt = T/n; }

    // Getter methods
    double getn() {return n;}
    double getdt() {return dt;}

    // Setter methods
    void setn(double n) { this->n = n; }
    void setdt(double dt) { this->dt = dt; }
};

class LookbackPut : public LookbackOption
{
public:
    // Constructor
    LookbackPut(double S, double K, double r, double sigma, double T, int n) :
        LookbackOption(S, K, r, sigma, T, n) {}

    // Function to price a lookback put option using Monte Carlo methods
    double price(int num_simulations) {
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
};


class LookbackCall : public LookbackOption
{
public:
    // Constructor
    LookbackCall(double S, double K, double r, double sigma, double T, int n) :
        LookbackOption(S, K, r, sigma, T, n) {}

    // Function to price a lookback call option using Monte Carlo methods
    double price(int num_simulations) {
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
};

























int main()
{
    double S = 100;      // Current price of the underlying asset
    double K = 110;      // Strike price
    double r = 0.05;     // Risk-free interest rate
    double sigma = 0.2;  // Volatility of the underlying asset
    double T = 1;        // Time to expiration (in years)
    int n = 12;
    int num_sim = 6000;

    // Create an European call option object
    EuropeanCall call(S, K, r, sigma, T);
    // Calculate the price of the call option
    double call_price = call.price();
    std::cout << "Price of European call option: " << call_price << std::endl;

    // Create an European put option object
    EuropeanPut put(S, K, r, sigma, T);
    // Calculate the price of the put option
    double put_price = put.price();
    std::cout << "Price of European put option: " << put_price << std::endl;

    call.replicate();

    double delta_call = call.delta();
    std::cout << "Delta of European Call option: " << delta_call << std::endl;

    AsianCall asiancall(S,K,r,sigma,T,n);
    double asian_call_price = asiancall.price(num_sim);
    std::cout << "Price of Asian call option: " << asian_call_price << std::endl;

    AsianPut asianput(S,K,r,sigma,T,n);
    double asian_put_price = asianput.price(num_sim);
    std::cout << "Price of Asian put option: " << asian_put_price << std::endl;


    LookbackCall lookbackcall(S,K,r,sigma,T,n);
    double lookbackcall_price = lookbackcall.price(num_sim);
    std::cout << "Price of Lookback call option: " << lookbackcall_price << std::endl;

    LookbackPut lookbackput(S,K,r,sigma,T,n);
    double lookbackput_price = lookbackput.price(num_sim);
    std::cout << "Price of Lookback put option: " << lookbackput_price << std::endl;



    return 0;
}



