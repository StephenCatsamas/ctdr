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
    out.log(INF) << "================ ctdr ===============" << std::endl;

    {  
    out.log(INF) << "loading phantom" << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);

    
    field<double> sinogram;    
    sinogram_intensity(phantom, projections, sinogram);
    // sinogram_attenuation(phantom, projections, sinogram);

    doub2png("data/sinogram.png", sinogram, SCALE_ZERO);
    sinogram.save("data/sinogram.fld");
    }
    {
    field<double> tomogram;
    field<double> sinogram = field<double>::load("data/sinogram.fld");

    intensity2attentuation(sinogram);
    doub2png("data/sinogram_att.png", sinogram, SCALE_ZERO);

    std::chrono::steady_clock::time_point tik,tok;
    
    tik = std::chrono::high_resolution_clock::now();
    recon_bp(sinogram, projections, tomogram);
    tok = std::chrono::high_resolution_clock::now();
    
    out.log(INF) << " bp: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    doub2png("data/tomogram_bp.png", tomogram);
    
    tik = std::chrono::high_resolution_clock::now();
    recon_fbp(sinogram, projections, tomogram);
    tok = std::chrono::high_resolution_clock::now();

    out.log(INF) << "fbp: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    doub2png("data/tomogram_fbp.png", tomogram);
    
    const int padding_factor = 8.0;
    tik = std::chrono::high_resolution_clock::now();
    recon_dfi(sinogram, projections, padding_factor, tomogram);
    tok = std::chrono::high_resolution_clock::now();
    
    out.log(INF) << "fdi: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    doub2png("data/tomogram_dfi.png", tomogram);
    
    tik = std::chrono::high_resolution_clock::now();
    recon_art(sinogram, projections, tomogram);
    tok = std::chrono::high_resolution_clock::now();
    
    out.log(INF) << "art: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    doub2png("data/tomogram_art.png", tomogram);
    }
    return 0;
}
