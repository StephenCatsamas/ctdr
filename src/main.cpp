#include <iostream>
#include <algorithm>
#include <string>
#include <complex>
#include <vector>

#include <chrono>

#include "field.h"
#include "util.h"
#include "recon.h"
#include "sim.h"

const int projections = 512;

int main() {
    out.log(INF) << "ctdr" << std::endl;
      
    out.log(INF) << "loading phantom" << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);
        
    // auto tomogram = field<double>(phantom.height, phantom.width);
    auto sinogram = field<double>(phantom.height, phantom.width);
    
    sinogram_intensity(phantom, projections, sinogram);

    doub2png("data/sinogram.png", sinogram);

    // std::chrono::steady_clock::time_point tik,tok;
    
    // tik = std::chrono::high_resolution_clock::now();
    // recon_bp(phantom, projections, tomogram);
    // tok = std::chrono::high_resolution_clock::now();
    
    // out.log(INF) << " bp: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    // doub2png("data/tomogram_bp.png", tomogram);
    
    // tik = std::chrono::high_resolution_clock::now();
    // recon_fbp(phantom, projections, tomogram);
    // tok = std::chrono::high_resolution_clock::now();

    // out.log(INF) << "fbp: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    // doub2png("data/tomogram_fbp.png", tomogram);
    
    // const int padding_factor = 8.0;
    // tik = std::chrono::high_resolution_clock::now();
    // recon_dfi(phantom, projections, padding_factor, tomogram);
    // tok = std::chrono::high_resolution_clock::now();
    
    // out.log(INF) << "fdi: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    // doub2png("data/tomogram_dfi.png", tomogram);
    
    // tik = std::chrono::high_resolution_clock::now();
    // recon_art(phantom, projections, tomogram);
    // tok = std::chrono::high_resolution_clock::now();
    
    // out.log(INF) << "art: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    // doub2png("data/tomogram_art.png", tomogram);
    
    return 0;
}
