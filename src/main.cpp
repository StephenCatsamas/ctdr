#include <iostream>

#include "field.h"

#include "util.h"

int project(const field<double>& phantom, const double angle, std::vector<double> projection){
    
    
    return 1;
}



int main() {
    out.log(INF) << "ctdr" << std::endl;
    
    out.log(INF) << "loading phantom" << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);
    
    phantom.rot(0.154);
    
    std::cout << "saving phantom" << std::endl;
    doub2png("data/output.png", phantom);
    
    return 0;
}
