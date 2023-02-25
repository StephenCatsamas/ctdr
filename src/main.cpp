#include <iostream>
#include <algorithm>
#include <string>
#include <complex>
#include <vector>

#include "field.h"
#include "util.h"
#include "recon.h"

const int projections = 256;

int main() {
    out.log(INF) << "ctdr" << std::endl;
      
    out.log(INF) << "loading phantom" << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);
        
    auto tomogram = field<double>(phantom.height, phantom.width);
        
    recon_bp(phantom, projections, tomogram);
    
    tomogram.clamp(0.0,1.0);
    out.log(INF) << "saving bp" << std::endl;
    doub2png("data/tomogram_bp.png", tomogram);
    
    recon_fbp(phantom, projections, tomogram);
    
    tomogram.clamp(0.0,1.0);
    out.log(INF) << "saving fbp" << std::endl;
    doub2png("data/tomogram_fbp.png", tomogram);
    
    recon_dfi(phantom, projections, tomogram);
    
    tomogram.clamp(0.0,1.0);
    out.log(INF) << "saving dfi" << std::endl;
    doub2png("data/tomogram_dfi.png", tomogram);
    
    recon_art(phantom, projections, tomogram);
    
    tomogram.clamp(0.0,1.0);
    out.log(INF) << "saving art" << std::endl;
    doub2png("data/tomogram_art.png", tomogram);
    
    
    

    return 0;
}
