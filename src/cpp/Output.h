/**
    Output header
    Output.h
    Purpose: Functions to print with colors and in logs.

    @author J.M.G. Salmerón
    @version 1.0 26/09/2018
    
    Created on: 26/09/2018
*/

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iostream>
#include <fstream>
#include <sstream>

enum Color {
    NONE = 0,
    BLACK, RED, GREEN,
    YELLOW, BLUE, MAGENTA,
    CYAN, WHITE
};

std::string color(Color foreground = NONE, Color background = NONE);

#endif /* OUTPUT_H_ */