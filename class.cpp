// Kai Rowlands
// Numerical Differentiation & Numerical Integration

#include <functional> 

using namespace std;

class FunctionTools{
public:
    std::function<double(double)> func; // Function that we will be applying the methods to

    FunctionTools(std::function<double(double)> x) : func(x){}  // Constructor: take the function and assign to variable func

    // Define the forward difference
    double ForwardDifference(double x, double h){
        return (func(x + h) - func(x)) / h;
    }

    // Define the backward difference
    double BackwardDifference(double x, double h){
        return (func(x) - func(x - h)) / h;
    }

    // Define the central difference
    double CentralDifference(double x, double h){
        return (func(x + h) - func(x - h)) / (2.0 * h);
    }

    // Define the left Riemann sum
    double LeftRiemann(double a, double b, int N){
        double xi;
        double h = (b - a) / N;                     // Define the length of the subintervals
        double I = 0;                               // Initialize the numerical integral to zero

        // Loop for each value of N
        for(int i = 0; i < N; i++){
            xi = a + i * h;         // Define the endpoints of each subinterval
            I = I + h * func(xi);   // Compute the area of the rectangle using the point x1
        }
        return I;   // Return the numerical integral
    }

    // Define the right Riemann sum
    double RightRiemann(double a, double b, int N){
        double xi;
        double h = (b - a) / N;                     // Define the length of the subintervals
        double I = 0;                               // Initialize the numerical integral to zero

        // Loop for each value of N
        for(int i = 0; i < N; i++){
            xi = a + (i + 1.0) * h; // Define the endpoints of each subinterval
            I = I + h * func(xi);   // Compute the area of the rectangle using the point x1
        }
        return I;   // Return the numerical integral
    }
    
    // Define the midpoint rule
    double midpoint(double a, double b, int N){
        double xi;
        double h = (b - a) / N;                 // Define the length of the subintervals
        double I = 0;                           // Initialize the numerical integral to zero

        // Loop for each value of N
        for(int i = 0; i < N; i++){
            xi = a + (i + 0.5) * h; // Define the midpoint of each subinterval
            I = I + h * func(xi);   // Compute the area of the rectangle using the point x1
        }
        return I;   // Return the numerical integral
    }
    
    // Define the trapezoidal rule
    double trapezoid(double a, double b, int N){
        double x1, x2;
        double h = (b - a) / N;                     // Define the length of the subintervals
        double I = 0;                               // Initialize the numerical integral to zero

        // Loop for each value of N
        for(int i = 0; i < N; i++){
            x1 = a + i * h;                             // Define the first point along the curve
            x2 = a + (i + 1.0) * h;                     // Define the second point along the curve
            I = I + h / 2.0 * (func(x1) + func(x2));    // Compute the area of the trapezoid using the points x1 and x2
        }
        return I;   // Return the numerical integral
    }

    // Define Simpson's rule
    double simpson(double a, double b, int N){
        double x1, x2, x3;
        double h = (b - a) / N;                 // Define the length of the subintervals
        double I = 0;                           // Initialize the numerical intergral to zero

        // Loop for each value of N
        for(int i = 0; i < N; i++){
            x1 = a + i * h;                                             // Define the first point along the curve
            x2 = a + (i + 0.5) * h;                                     // Define the second point along the curve (the midpoint of x1 and x3)
            x3 = a + (i + 1.0) * h;                                     // Define the third point along the curve
            I = I + h / 6.0 * (func(x1) + 4.0 * func(x2) + func(x3));   // Compute Simpson's rule using the points x1, x2, and x3
        }
        return I;   // Return the numerical integral
    }
};