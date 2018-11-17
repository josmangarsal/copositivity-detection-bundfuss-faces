/**
    Output header
    Output.h
    Purpose: Functions to print with colors and in logs.

    @author J.M.G. Salmerón
    @version 1.0 26/09/2018
    
    Created on: 26/09/2018
*/

#include "Output.h"

std::string color(Color foreground, Color background) {
    std::stringstream s;
    s << "\033[";
    if (!foreground && ! background) {
        s << "0"; // reset colors if no params
    }
    if (foreground) {
        s << 29 + foreground;
        if (background) s << ";";
    }
    if (background) {
        s << 39 + background;
    }
    s << "m";
    return s.str();
}