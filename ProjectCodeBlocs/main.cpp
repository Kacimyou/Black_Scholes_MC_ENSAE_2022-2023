/*
	Paul Wattellier, Mira Maamari, Denisa Draghian, Kacim Younsi
	27/01/2022
	Assignment C++ - End of term Project
	Main
*/

#include <stdio.h>
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
    double T = 1;  // Time to expiration (in years)
    double T1=0.5;// Time to exercise the chooser option
    int n = 12;
    int num_sim = 6000;
    double K1= 113; // the second strike price for Gap Options
    double barrier=98; //the barrier for Barrier Option
    int type=1; // the type of barrier call: in or out
    int type1=0; // the type of barrier call: up or down


    // Create an European call option object
    EuropeanCall call(S, K, r, sigma, T);
    // Calculate the price of the call option
    double call_price = call.price();
    std::cout << "Price of European call option: " << call_price << std::endl;

    double call_price_1 = call.price_MonteCarlo(6000);
    std::cout << "Price of European call option using Monte Carlo simulation: " << call_price_1 << std::endl;

    double difference=call.difference(6000);
    std::cout << "Difference of prices: " << difference<<std::endl;



    // Create an European put option object
    EuropeanPut put(S, K, r, sigma, T);
    // Calculate the price of the put option
    double put_price = put.price();
    std::cout << "Price of European put option: " << put_price << std::endl;

    double put_price_1 = put.price_MonteCarlo(10000);
    std::cout << "Price of European put option using Monte Carlo simulation: " << put_price_1 << std::endl;


    call.replicate();
    put.replicate();

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

    BinaryCall binarycall(S,K,r,sigma,T,n);
    double binary_call_price = binarycall.price(num_sim);
    std::cout << "Price of Binary call option: " << binary_call_price << std::endl;

    BinaryPut binaryput(S,K,r,sigma,T,n);
    double binary_put_price = binaryput.price(num_sim);
    std::cout << "Price of Binary put option: " << binary_put_price << std::endl;

    ChooserOption chooser(S,K,r,sigma,T,T1);
    double chooser_price = chooser.price();
    std::cout << "Price of Chooser option: " << chooser_price << std::endl;

    GapCall gapcall(S,K,r,sigma,T,K1);
    double gap_call_price = gapcall.price();
    std::cout << "Price of Gap call option: " << gap_call_price << std::endl;

    BarrierCall barriercall(S,K,r,sigma,T,n,type,barrier,type1);
    double barrier_call_price= barriercall.price(num_sim);
    std::cout << "Price of Barrier call option: " << barrier_call_price << std::endl;


    return 0;
}
