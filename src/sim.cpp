
#include "sim.h"

#include <random>

enum integrand{
    INTENSITY,
    ATTENUATION,
};

int project_stem(const field<double>& phantom, const double angle, std::vector<double>& projection, integrand i_mode){
    
    int seed = 42;
    std::default_random_engine generator(seed);
    std::poisson_distribution<int> poission(1);
    int I0 = 1E3; //photons
    double proj_att_sf = 1.0/phantom.height;//attenuation factor to make sure that projections dont get toooo small in the exponent we will limit to 1/e.

    auto oriented = phantom;
    
    oriented.rot(angle);
    
    projection.resize(oriented.width);
    std::fill(projection.begin(), projection.end(), 0.0);
    
    for (int j = 0; j< oriented.width; j++){
        for (int i = 0; i< oriented.height; i++){
            projection[j] += oriented[i][j];
        }
        switch(i_mode){
            case INTENSITY:
                projection[j] = std::round(I0*exp(-projection[j]*proj_att_sf));
                break;
            case ATTENUATION:
                break;
        }
    }
    
    
    return 1;
}

int project_intensity(const field<double>& phantom, const double angle, std::vector<double>& projection){
    return project_stem(phantom, angle, projection, INTENSITY);
}

int project_attenuation(const field<double>& phantom, const double angle, std::vector<double>& projection){
    return project_stem(phantom, angle, projection, ATTENUATION);
}

int sinogram_base(const field<double>& phantom, int projections, field<double>& sinogram, integrand i_mode){
    double angle_step = 2*pi/projections;
    
    sinogram = field<double>(projections, phantom.width, 0.0);

    auto projection = std::vector<double>(phantom.width);
    for(int i = 0; i < projections; i++){
        double angle = i*angle_step;
        project_stem(phantom, angle, projection, i_mode);
        
        std::copy(projection.begin(), projection.end(), 
                  sinogram.data.begin()+(sinogram.width*i));
        
    }
    return 1;
}

int sinogram_intensity(const field<double>& phantom, int projections,field<double>& sinogram){
    return sinogram_base(phantom, projections, sinogram, INTENSITY);
}

int sinogram_attenuation(const field<double>& phantom, int projections,field<double>& sinogram){
    return sinogram_base(phantom, projections, sinogram, ATTENUATION);
}