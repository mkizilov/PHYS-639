#include <iostream>
#include <limits>
#include <cmath>
#include <random>
#include <iomanip>
#include "TH1F.h"
#include "TCanvas.h"

// #include <boost/math/distributions/gamma.hpp>

// Declare functions
long double Functions(int func, long double x);
std::string FunctionNames(int func);
bool equal(long double x, long double y);
long double MidpointIntegral(int func, long double a, long double b, int bins);
long double TrapezoidalIntegral(int func, long double a, long double b, int bins);
long double MonteCarloIntegral(int func, long double a, long double b, int bins);
long double MonteCarloRootMeanSquare(long double a, long double b, int bins, int func);
long double MonteCarloUncertainty(long double a, long double b, int bins, int func);
void plotFunctions(int func, long double a, long double b, int maxBins, int binRange);
void plotRootMeanSquare(int func, long double a, long double b, int maxBins, int binRange);
void plotUncertainty(int func, long double a, long double b, int maxBins, int binRange);


long double Functions(int func, long double x)
{
    if (func == 0) return 1; // for function f = 1
    else if (func == 1) return x; // for function f = x
    else if (func == 2) return (x*x*x); // for function f = x^2
    else if (func == 3) return (x*x)/((1+x*x)*(1+x*x)); // for function f = x^2/(1+x^2)^2
    else if (func == 4) return (300*exp(-pow((x-2),2)/0.00001)+3*exp(x-2)); // for function f = 300*exp(-pow((x-2),2)/0.00001)+3*exp(x-2)
    else if (func == 5) return (sin(x)*sin(x));
    else return 0; // for invalid function
}

std::string FunctionNames(int func)
{
    if (func == 0) return "f = 1";
    else if (func == 1) return "f = x";
    else if (func == 2) return "f = x^3";
    else if (func == 3) return "f = x^2/(1+x^2)^2";
    else if (func == 4) return "f = 300*exp(-pow((x-2),2)/0.00001)+3*exp(x-2)";
    else if (func == 5) return "f = sin(x)*sin(x)";
    else return "Invalid function";
}

// A function to check if two numbers are equal within the tolerance
bool equal(long double x, long double y)
{
    return std::abs(x - y) <= std::numeric_limits<long double>::epsilon();
}

long double MonteCarloRootMeanSquare(long double a, long double b, int bins, int func)
{
    long double sum = 0;
    long double integral = MonteCarloIntegral(func, a, b, bins);
    for (int i = 0; i < bins; i++)
    {
        long double x = a + (b-a)*i/bins;
        sum += pow(Functions(func, x) - integral, 2);
    }
    return sqrt(sum/bins);
}
long double MonteCarloUncertainty(long double a, long double b, int bins, int func)
{
    return MonteCarloRootMeanSquare(a, b, bins, func)/sqrt(bins);
}

long double MidpointIntegral(int func, long double a, long double b, int bins)
{
    // Check if the lower limit is infinity
    if (std::isinf(a) && a < 0)
    {
        // Use a large negative number instead
        a = -1e9;
    }
    // Check if the upper limit is infinity
    if (std::isinf(b) && b > 0)
    {
        // Use a large positive number instead
        b = 1e9;
    }
    long double sum = 0;
    long double dx = (b-a)/bins;
    for (int i = 0; i < bins; i++)
    {
        long double middle = a + (i*dx) + (dx*0.5); // Midpoint of each bin
        sum += Functions(func, middle)*dx;
    }
    return sum;
}

long double TrapezoidalIntegral(int func, long double a, long double b, int bins)
{
    // Check if the lower limit is infinity
    if (std::isinf(a) && a < 0)
    {
        // Use a large negative number instead
        a = -1e9;
    }
    // Check if the upper limit is infinity
    if (std::isinf(b) && b > 0)
    {
        // Use a large positive number instead
        b = 1e9;
    }
    long double sum = 0;
    long double dx = (b-a)/bins;
    for (int i = 0; i < bins; i++)
    {
        long double x0 = a + (i*dx);
        long double x1 = a + ((i+1)*dx);
        sum += (Functions(func, x0) + Functions(func, x1))*dx*0.5;
    }
    return sum;
}





long double MonteCarloIntegral(int func, long double a, long double b, int bins)
{
    long double sum1 = 0;
    long double sum2 = 0;
    long double sum = 0;
    long double x1, x2;
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<long double> dis(a, b);
    long double x;
    std::normal_distribution<long double> dis1(2, 1);
    std::gamma_distribution<long double> dis2(2, 1);

// Divide the function value by the sampling density
    
    // Generate random number from the chosen PDF
    for (int i = 0; i < bins; i++)
    {   
        x = dis(gen);
        x1 = dis1(gen);
        x2 = dis2(gen);
        while (x1 < a || x1 > b)
        {
            x1 = dis1(gen);
        }
        while (x2 < a || x2 > b)
        {
            x2 = dis2(gen);
        }
        //std::cout << "x = " << x1 << " f(x) = " << Functions(func, x1) << " pdf = " << ROOT::Math::normal_pdf(x1, 0.5, 2.5) << std::endl;
        sum1 += Functions(func, x1) ;
        sum2 += Functions(func, x2) / ROOT::Math::gamma_pdf(x2, 1, 2);
        sum += Functions(func, x);
        
        //std::cout << "x = " << x << " f(x) = " << Functions(func, x) << " pdf = " << ROOT::Math::normal_pdf(x, 0.1, 2) << std::endl;
    }

    std::cout << "sum1 = " << sum1*(b-a)/bins << " sum2 = " << sum2*(b-a)/bins << std::endl;
    std::cout << "sum = " << sum*(b-a)/bins << std::endl;
    return sum/bins*(b-a);
    
}








void plotFunctions(int func, long double a, long double b, int maxBins, int binRange)
{
    // Create a new canvas
    TCanvas *c1 = new TCanvas("c1", "Integral Convergence", 800, 600);

    // Create TGraphs for each method
    TGraph *gMidpoint = new TGraph();
    TGraph *gTrapezoidal = new TGraph();
    TGraph *gMonteCarlo = new TGraph();

    // Set line colors for each graph
    gMidpoint->SetLineColor(kRed);
    gTrapezoidal->SetLineColor(kBlue);
    gMonteCarlo->SetLineColor(kGreen);
    // Calculate the integral for each method and add the points to the graphs
    int pointIndex = 0;
    for (int bins = 1; bins <= maxBins; bins = bins + binRange)
    {
        long double MIintegral = MidpointIntegral(func, a, b, bins);
        long double Tintegral = TrapezoidalIntegral(func, a, b, bins);
        long double MCintegral = MonteCarloIntegral(func, a, b, bins);

        gMidpoint->SetPoint(pointIndex, bins, MIintegral);
        gTrapezoidal->SetPoint(pointIndex, bins, Tintegral);
        gMonteCarlo->SetPoint(pointIndex, bins, MCintegral);

        pointIndex++;
    }

    // Set marker style for each graph
    gMidpoint->SetMarkerStyle(20);
    gTrapezoidal->SetMarkerStyle(20);
    gMonteCarlo->SetMarkerStyle(20);

    // Set X log scale
    gPad->SetLogx();


    // Draw the graphs on the canvas
    gMidpoint->Draw("AL");
    gTrapezoidal->Draw("same");
    gMonteCarlo->Draw("same");

    // X log scale


    // Add a legend top right
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(gMidpoint,"Midpoint","l");
    legend->AddEntry(gTrapezoidal,"Trapezoidal","l");
    legend->AddEntry(gMonteCarlo,"Monte Carlo","l");
    legend->Draw();

    // Calculate convergence rate
    long double MIintegral = MidpointIntegral(func, a, b, maxBins);
    long double Tintegral = TrapezoidalIntegral(func, a, b, maxBins);
    long double MCintegral = MonteCarloIntegral(func, a, b, maxBins);
    std::cout << "The Midpoint Integral of the function " << FunctionNames(func) <<" from " << a << " to " << b << " is " << MIintegral << std::endl;
    std::cout << "The Trapezoidal Integral of the function " << FunctionNames(func) <<" from " << a << " to " << b << " is " << Tintegral  << std::endl;
    std::cout << "The Monte Carlo Integral of the function " << FunctionNames(func) <<" from " << a << " to " << b << " is " << MCintegral << std::endl;
    std::cout << "The convergence rate of the Midpoint Integral of the function " << FunctionNames(func) <<" from " << a << " to " << b << " is " << MIintegral - MidpointIntegral(func, a, b, maxBins/2) << std::endl;
    std::cout << "The convergence rate of the Trapezoidal Integral of the function " << FunctionNames(func) <<" from " << a << " to " << b << " is " << Tintegral - TrapezoidalIntegral(func, a, b, maxBins/2) << std::endl;
    std::cout << "The convergence rate of the Monte Carlo Integral of the function " << FunctionNames(func) <<" from " << a << " to " << b << " is " << MCintegral - MonteCarloIntegral(func, a, b, maxBins/2) << std::endl;

    // Set y range to be the same for all graphs
    gMidpoint->GetYaxis()->SetRangeUser(0.9*gMidpoint->GetHistogram()->GetMinimum(), 1.1*gMidpoint->GetHistogram()->GetMaximum());
    string title = "Convergence of " + (FunctionNames(func)) + " from " + std::to_string(a) + " to " + std::to_string(b);
    // Add titles and labels
    gMidpoint->SetTitle(title.c_str());
    gMidpoint->GetXaxis()->SetTitle("Bins or Iterations");
    gMidpoint->GetYaxis()->SetTitle("Integral Value");



    // Update the canvas
    c1->Update();
    string filename = (to_string(func))+"_simple_convergence.pdf";
    // Save the plot
    c1->SaveAs(filename.c_str());
}

void plotRootMeanSquare(int func, long double a, long double b, int maxBins, int binRange)
{
    // Create a new canvas
    TCanvas *c1 = new TCanvas("c1", "Root Mean Square Convergence", 800, 600);

    // Create TGraphs for each method
    TGraph *gMidpoint = new TGraph();

    // Set line colors for each graph
    gMidpoint->SetLineColor(kRed);

    // Calculate the integral for each method and add the points to the graphs
    int pointIndex = 0;
    for (int bins = 1; bins <= maxBins; bins = bins + binRange)
    {
        long double MIintegral = MonteCarloRootMeanSquare(a, b, bins, func);
        gMidpoint->SetPoint(pointIndex, bins, MIintegral);
        pointIndex++;
    }

    // Set marker style for each graph
    gMidpoint->SetMarkerStyle(20);

    // Draw the graphs on the canvas
    gMidpoint->Draw("AL");

    // Add a legend top right
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetFillStyle(0);
    legend->AddEntry(gMidpoint,"Root Mean Square Error MC","l");
    legend->Draw();

    // Set log scale
    gPad->SetLogx();

    // Add titles and labels
    string title = "Root Mean Square Error Convergence of " + (FunctionNames(func)) + " from " + std::to_string(a) + " to " + std::to_string(b);
    gMidpoint->SetTitle(title.c_str());
    gMidpoint->GetXaxis()->SetTitle("Bins or Iterations");
    gMidpoint->GetYaxis()->SetTitle("RMSE");

    // Update the canvas
    c1->Update();

    // Save the plot
    string filename = (to_string(func))+"_RMS_convergence.pdf";
    c1->SaveAs(filename.c_str());

}

void plotUncertainty(int func, long double a, long double b, int maxBins, int binRange)
{
    // Create a new canvas
    TCanvas *c1 = new TCanvas("c1", "Uncertainty Convergence", 800, 600);

    // Create TGraphs for each method
    TGraph *gMidpoint = new TGraph();

    // Set line colors for each graph
    gMidpoint->SetLineColor(kRed);

    // Calculate the integral for each method and add the points to the graphs
    int pointIndex = 0;
    for (int bins = 1; bins <= maxBins; bins = bins + binRange)
    {
        long double MIintegral = MonteCarloUncertainty(a, b, bins, func);
        gMidpoint->SetPoint(pointIndex, bins, MIintegral);
        pointIndex++;
    }

    // Set marker style for each graph
    gMidpoint->SetMarkerStyle(20);

    // Draw the graphs on the canvas
    gMidpoint->Draw("AL");

    // Add a legend top right
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(gMidpoint,"Uncertainty MC","l");
    legend->Draw();

    // Set log scale
    gPad->SetLogx();

    // Add titles and labels
    string title = "Uncertainty Convergence of " + (FunctionNames(func)) + " from " + std::to_string(a) + " to " + std::to_string(b);
    gMidpoint->SetTitle(title.c_str());
    gMidpoint->GetXaxis()->SetTitle("Bins or Iterations");
    gMidpoint->GetYaxis()->SetTitle("Uncertainty");

    // Update the canvas
    c1->Update();
    string filename = (to_string(func))+"_uncertainty_convergence.pdf";
    // Save the plot
    c1->SaveAs(filename.c_str());

}


void MK_LAB01(int func = 2, std::string lowerLimit = "1", std::string upperLimit = "3", int bins = 100)
{
    long double a = std::stold(lowerLimit);
    long double b;
    if (upperLimit == "pi2")
    {
     b = M_PI/2;
    }
    else b = std::stold(upperLimit);
    long double MIintegral = MidpointIntegral(func, a, b, bins);
    long double Tintegral = TrapezoidalIntegral(func, a, b, bins);
    long double MCintegral = MonteCarloIntegral(func, a, b, bins);

    std::cout << "The Midpoint Integral of the function " << FunctionNames(func) <<" from " << a << " to " << b << " is " << MIintegral << std::endl;
    std::cout << "The Trapezoidal Integral of the function " << FunctionNames(func) <<" from " << a << " to " << b << " is " << Tintegral  << std::endl;
    std::cout << "The Monte Carlo Integral of the function " << FunctionNames(func) <<" from " << a << " to " << b << " is " << MCintegral << std::endl;
    
    // plotFunctions(func, a, b, bins, 1);
    // plotRootMeanSquare(func, a, b, bins, 1);
    // plotUncertainty(func, a, b, bins, 1);
}
