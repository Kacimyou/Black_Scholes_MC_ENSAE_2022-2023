/*
	Paul Wattellier, Mira Maamari, Denisa Draghian, Kacim Younsi
	27/01/2022
	Assignment C++ - End of term Project
	Main
*/


#include "OptionPricing.h"
#include <cmath>
#include <math.h>
#include <iostream>
#include <random>



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
