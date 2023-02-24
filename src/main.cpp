#include <iostream>
#include <numbers>
#include <algorithm>
#include <string>
#include <complex>
#include <vector>

#include "field.h"

#include "util.h"

#include "pocketfft.h"

using namespace pocketfft;

double pi = std::numbers::pi_v<double>;

#define PROJECTIONS (256)

int project(const field<double>& phantom, const double angle, std::vector<double>& projection){
    
    auto oriented = phantom;
    
    oriented.rot(angle);
    
    projection.resize(oriented.width);
    std::fill(projection.begin(), projection.end(), 0.0);
    
    for (int j = 0; j< oriented.width; j++){
    for (int i = 0; i< oriented.height; i++){
        projection[j] += oriented[i][j];
    }
    }
    
    return 1;
}

int proj2png(char* fp, std::vector<double>& projection){
    const int projection_plot_height = 256;
    
    const int projection_size = projection.size();
    
    auto plot = field<double>(projection_plot_height, projection_size, 1);
    
    double mx = *std::max_element(projection.begin(), projection.end());
    
    for(int j = 0; j < projection_size; j++){
        int i = floor((projection_plot_height-1)*(1.0-projection[j]/mx));
        plot[i][j] = 0;
    }
    doub2png(fp, plot);
    return 1;
}

int back_project(const std::vector<double>& projection, double angle, field<double>& back_projection){
    
    for (int j = 0; j< back_projection.width; j++){
    for (int i = 0; i< back_projection.height; i++){
        back_projection[i][j] = projection[j];
    }
    } 
    back_projection.rot(-angle);
    
    return 1;
}

int recon_bp(const field<double>& phantom, field<double>& tomogram){
    
    const double angle_step = 2.0*pi/PROJECTIONS;
    
    tomogram.fill(0.0);
    for(double angle = 0.0; angle < 2*pi; angle += angle_step){
        out.log(INF) << angle << std::endl;
        auto back_projection = field<double>(tomogram.height,tomogram.width);
        std::vector<double> projection;
        project(phantom, angle, projection);
        back_project(projection, angle, back_projection);
        tomogram += back_projection* (0.5/(PROJECTIONS*back_projection.height));
       
    }
    return 1; 
}

int recon_fbp(const field<double>& phantom, field<double>& tomogram){
    
    const double angle_step = 2.0*pi/PROJECTIONS;   
    
    tomogram.fill(0.0);
    for(double angle = 0.0; angle < 2.0*pi; angle += angle_step){
        out.log(INF) << angle << std::endl;
        auto back_projection = field<double>(tomogram.height,tomogram.width);
        std::vector<double> projection;
        project(phantom, angle, projection);
        
        //fft stuff
        //pad projection
        projection.resize(2*projection.size(), 0.0);

        shape_t shape = {projection.size()};
        stride_t stride_r = {sizeof(double)};
        stride_t stride_c = {sizeof(std::complex<double>)};
        shape_t axes = {0};

        auto f_projection = std::vector<std::complex<double>>(shape[0]/2 + 1);
        double scale_factor = 0.5/projection.size();

        r2c(shape,
            stride_r, 
            stride_c, 
            axes,
            FORWARD, 
            projection.data(), 
            f_projection.data(), 
            scale_factor);
        
        //filter
        for(int i = 0; i < shape[0]/2 + 1; i++){
            double weight = 2*pi*(i);
            // double window = 0.54+0.46*cos(2*pi*i/(4*shape[0]/2));
            f_projection[i] *= weight;
        }    
        
        c2r(shape,
            stride_c, 
            stride_r, 
            axes,
            BACKWARD, 
            f_projection.data(), 
            projection.data(), 
            1.0);
   
        //unpad projection
        projection.resize(back_projection.width);
        //end fft stuff    
                
        back_project(projection, angle, back_projection);
        tomogram += back_projection * (0.5/(PROJECTIONS*back_projection.height));
       
    }
    return 1; 
}

int main() {
    out.log(INF) << "ctdr" << std::endl;
    
    out.log(INF) << "loading phantom" << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);
    
    auto tomogram = field<double>(phantom.height, phantom.width);
        
    // recon_bp(phantom, tomogram);
    recon_fbp(phantom, tomogram);
    
    double tm_mx = tomogram.max();
    out.log(INF) << tm_mx << std::endl;
    tomogram.clamp(0.0,1.0);
    
    out.log(INF) << "saving tomogram" << std::endl;
    doub2png("data/tomogram.png", tomogram);

    return 0;
}
