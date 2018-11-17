/**
    Compare body
    Compare.h
    Purpose: Functions to compare according to an accuracy.

    @author J.M.G. Salmerón
    @version 1.0 20/09/2018
    
    Created on: 20/09/2018
*/

#include "Compare.h"

bool eq(double a, double b) {
	if (a <= b && b - a < accuracy)
		return true;

	if (a >= b && a - b < accuracy)
		return true;

	return false;
}

bool lt(double a, double b) {
    if (eq(a, b))
		return false;

	if (a < b)
		return true;

	return false;
}

bool gt(double a, double b) {
    if (eq(a, b))
		return false;

	if (a > b)
		return true;

	return false;
}

bool le(double a, double b) {
    if (eq(a, b) || a < b)
		return true;

	return false;
}

bool ge(double a, double b) {
    if (eq(a, b) || a > b)
		return true;

	return false;
}