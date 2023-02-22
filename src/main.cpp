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


const int projection_plot_height = 256;
int main() {
    out.log(INF) << "ctdr" << std::endl;
    
    out.log(INF) << "loading phantom" << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);
    
    std::vector<double> projection;
    
    for(double angle = 0.0; angle < pi; angle += pi/12){
    // double angle = 0.1;
    
    project(phantom, angle, projection);
    
    auto plot = field<double>(projection_plot_height, phantom.width, 1);
    
    double mx = *std::max_element(projection.begin(), projection.end());
    
    for(int j = 0; j < phantom.width; j++){
        int i = floor((projection_plot_height-1)*(1.0-projection[j]/mx));
        plot[i][j] = 0;
    }
    
    char proj_fp[128];
    sprintf(proj_fp, "data/proj%f.png", angle);
    
    std::cout << "saving phantom" << std::endl;
    doub2png("data/output.png", phantom);
    doub2png("data/projection.png", plot);
    doub2png(proj_fp, plot);
    }
    return 0;
}
