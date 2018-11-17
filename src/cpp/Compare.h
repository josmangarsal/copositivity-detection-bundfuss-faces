/**
    Compare header
    Compare.h
    Purpose: Functions to compare according to an accuracy.

    @author J.M.G. Salmerón
    @version 1.0 20/09/2018
    
    Created on: 20/09/2018
*/

#ifndef COMPARE_H_
#define COMPARE_H_

#define accuracy 10e-6
#define Chol_accuracy 10e-6

bool eq(double a, double b); // a == b
bool lt(double a, double b); // a < b
bool gt(double a, double b); // a > b
bool le(double a, double b); // a <= b
bool ge(double a, double b); // a >= b

#endif /* COMPARE_H_ */
