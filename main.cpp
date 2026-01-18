// Kai Rowlands
// Numerical Differentiation & Numerical Integration

#include "class.cpp"
#include <cmath>        // Used for the absolute value function and exponentiation
#include <iomanip>      // Used for scientific, setprecision, and setw
#include <iostream>

using namespace std;

// Define the cubic function
double cubic(double x){
    return pow(x, 3.0);
}

// Define the derivative of cubic
double dcubic(double x){
    return 3.0 * pow(x, 2.0);
}

// Define the antiderivative of cubic
double acubic(double a, double b){
    return (pow(b, 4.0) - pow(a, 4.0)) / 4.0;
}

// Define the trigonometric function 
double trig(double x){
    return sin(5.0 * x) * cos(2.0 * x);
}

// Define the derivative of trig
double dtrig(double x){
    return 5.0 * cos(5.0 * x) * cos(2.0 * x) - 2.0 * sin(5.0 * x) * sin(2.0 * x);
}

// Define the antiderivative of trig
double atrig(double a, double b){
    return -1.0 / 14.0 * cos(7.0 * b) - 1.0 / 6.0 * cos(3.0 * b) + 1.0 / 14.0 * cos(7.0 * a) + 1.0 / 6.0 * cos(3.0 * a);
}

int main(){
    double h;
    double errorF, errorB, errorC;                  // Errors of the forward difference, the backward difference, and the central difference
    double errorL, errorR, errorM, errorT, errorS;  // Errors of the left Riemann sum, the right Riemann sum, the midpoint rule, the trapezoidal rule, and Simpson's rule
    double df, af;                                  // Derivate of f and antiderivative of f

    FunctionTools f1(cubic);    // Objects of the class FunctionTools
    FunctionTools f2(trig);     // Change f2 to f1 below to perform the numerical methods on f1
    
    double x = 0.123;   // x-value where we are approximating the nuerical derivatives at
    
    cout << scientific << setprecision(3);                                                                          // Use scientific notation with three digits of accuracy
    cout << setw(12) << "h" << setw(12) << "Forward" << setw(12) << "Backward" << setw(12) << "Centered" << endl;   // Titles for each column
    
    // Loop for each value of h
    for(int i = 1; i < 8; i++){
        h = 1.0 / pow(2.0, i);    // Spacing between the points for the numerical derivative
        
        df = f2.ForwardDifference(x,h);     // Compute the forward difference
        errorF = abs(df - dtrig(x));        // Compute the absolute error of the forward difference
        df = f2.BackwardDifference(x,h);    // Compute the backward difference
        errorB = abs(df - dtrig(x));        // Compute the absolute error of the backward difference
        df = f2.CentralDifference(x,h);     // Compute the centeral difference
        errorC = abs(df - dtrig(x));        // Compute the absolute error of the central difference
        
        cout << setw(12) << h << setw(12) << errorF << setw(12) << errorB << setw(12) << errorC << endl;   // Reports the errors of the three numerical derivatives
    }
    cout << endl << endl;                                                                                                                                         // Creates an empty line to split up the charts between the numerical derivatives and the numerical integration methods
    cout << setw(12) << "N" << setw(12) << "Left" << setw(12) << "Right" << setw(12) << "Midpoint" << setw(12) << "Trapezoid" << setw(12) << "Simpson" << endl;   // Titles for each column
    
    // Loop for each value of N
    for(int i = 1; i < 8; i++){
        int N = pow(2.0, i);    // Number of subintervals for each numerical integration method
        
        af = f2.LeftRiemann(-1.0, 3.0, N);    // Compute the left Riemann sum
        errorL = abs(af - atrig(-1.0, 3.0));  // Compute the absolute error of the left Riemann sum
        af = f2.RightRiemann(-1.0, 3.0, N);   // Compute the right Riemann sum
        errorR = abs(af - atrig(-1.0, 3.0));  // Compute the absolute error of the right Riemann sum
        af = f2.midpoint(-1.0, 3.0, N);       // Compute the midpoint rule
        errorM = abs(af - atrig(-1.0, 3.0));  // Compute the absolute error of midpoint rule
        af = f2.trapezoid(-1.0, 3.0, N);      // Compute the trapezoid rule
        errorT = abs(af - atrig(-1.0, 3.0));  // Compute the absolute error of the trapezoidal rule
        af = f2.simpson(-1.0, 3.0, N);        // Compute Simpson's rule
        errorS = abs(af - atrig(-1.0, 3.0));  // Compute the absolute error of Simpson's rule
        
        cout << setw(12) << N << setw(12) << errorL << setw(12) << errorR << setw(12) << errorM << setw(12) << errorT << setw(12) << errorS << endl;    // Reports the errors of the five numerical integration methods
    }
    return 0;
}