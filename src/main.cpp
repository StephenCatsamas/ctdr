#include <iostream>
#include <numbers>
#include <algorithm>
#include <string>

#include "field.h"

#include "util.h"

double pi = std::numbers::pi_v<double>;

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
    
    const double angle_step = 2*pi/32;
    
    tomogram.fill(0.0);
    for(double angle = 0.0; angle < 2*pi; angle += angle_step){
        out.log(INF) << angle << std::endl;
        auto back_projection = field<double>(tomogram.height,tomogram.width);
        std::vector<double> projection;
        project(phantom, angle, projection);
        back_project(projection, angle, back_projection);
        tomogram += back_projection;
       
    }
    return 1; 
}

int main() {
    out.log(INF) << "ctdr" << std::endl;
    
    out.log(INF) << "loading phantom" << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);
    
    auto tomogram = field<double>(phantom.height, phantom.width);
        
    recon_bp(phantom, tomogram);
    
    double tm_mx = tomogram.max();
    tomogram *= 1.0/tm_mx;
       
    
    out.log(INF) << "saving tomogram" << std::endl;
    doub2png("data/tomogram.png", tomogram);

    return 0;
}
