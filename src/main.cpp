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

 scan_prop scan = {
        INTENSITY,//projection mode
        42,       //random seed
        10'000,   //I0
        512,      //projections
        1.0/256.0,//att_sf
        true,//quanisation noise
        true,
    };

std::chrono::steady_clock::time_point tik,tok;

int main() {
    out.log(INF) << "================ ctdr ===============" << std::endl;

    {  
    out.log(INF) << "SIMULATING PHANTOM..." << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);

    field<double> sinogram;  

    tik = std::chrono::high_resolution_clock::now();
    sinogram_sim(phantom, scan, sinogram);
    tok = std::chrono::high_resolution_clock::now();

    out.log(INF) << " sim: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    doub2png("data/sinogram.png", sinogram, SCALE_ZERO);
    sinogram.save("data/sinogram.fld");



    }

    {
    out.log(INF) << "RECONSTRUCTING PHANTOM..." << std::endl;
    field<double> tomogram;
    field<double> sinogram = field<double>::load("data/sinogram.fld");

    intensity2attentuation(sinogram, scan);

    
    
    tik = std::chrono::high_resolution_clock::now();
    recon_bp(sinogram, scan, tomogram);
    tok = std::chrono::high_resolution_clock::now();
    
    out.log(INF) << " bp: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    doub2png("data/tomogram_bp.png", tomogram);
    
    tik = std::chrono::high_resolution_clock::now();
    recon_fbp(sinogram, scan, tomogram);
    tok = std::chrono::high_resolution_clock::now();

    out.log(INF) << "fbp: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    doub2png("data/tomogram_fbp.png", tomogram);
    
    const int padding_factor = 8.0;
    tik = std::chrono::high_resolution_clock::now();
    recon_dfi(sinogram, scan, padding_factor, tomogram);
    tok = std::chrono::high_resolution_clock::now();
    
    out.log(INF) << "fdi: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    doub2png("data/tomogram_dfi.png", tomogram);
    
    tik = std::chrono::high_resolution_clock::now();
    recon_art(sinogram, scan, tomogram);
    tok = std::chrono::high_resolution_clock::now();
    
    out.log(INF) << "art: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik) << std::endl;
    doub2png("data/tomogram_art.png", tomogram);
    }
    return 0;
}
