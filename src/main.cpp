#include <iostream>

#include "field.h"

#include "util.h"


int main() {
    out.log(INF) << "ctdr" << std::endl;
    
    out.log(INF) << "loading phantom" << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);
    
    phantom[128][128] = 500;
    
    std::cout << "saving phantom" << std::endl;
    doub2png("data/output.png", phantom);
    
    return 0;
}
