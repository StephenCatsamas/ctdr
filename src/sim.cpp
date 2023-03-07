
#include "sim.h"

enum integrand{
    INTENSITY,
    ATTENUATION,
};

int project_stem(const field<double>& phantom, const double angle, std::vector<double>& projection, integrand i_mode){
    
    auto oriented = phantom;
    
    oriented.rot(angle);
    
    projection.resize(oriented.width);
    std::fill(projection.begin(), projection.end(), 0.0);
    
    for (int j = 0; j< oriented.width; j++){
    for (int i = 0; i< oriented.height; i++){
        swith(i_mode){
            case INTENSITY:
                projection[j] += exp(-oriented[i][j]);
                break;
            case ATTENUATION:
                projection[j] += oriented[i][j];
                break;
            
        }
    }
    }
    
    return 1;
}

int project_intensity(const field<double>& phantom, const double angle, std::vector<double>& projection){
    project_stem(phantom, angle, projection, INTENSITY);
}

int project_attenuation(const field<double>& phantom, const double angle, std::vector<double>& projection){
    project_stem(phantom, angle, projection, INTENSITY);
}

int sinogram_base(const field<double>& phantom, int projections, field<double>& sinogram, integrand i_mode){
    double angle_step = 2*pi/projections;
    
    sinogram = field<double>(projections, phantom.width, 0.0);
    auto projection = std::vector<double>(phantom.width);
    for(int i = 0; i < projections; i++){
        double angle = i*angle_step;
        project_intensity(phantom, angle, projection);
        
        std::copy(projection.begin(), projection.end(), 
                  sinogram.data.begin()+(sinogram.width*i));
        
    }
    return 1;
}

int sinogram_intensity(const field<double>& phantom, int projections,field<double>& sinogram){
    sinogram_base(phatom, projections, sinogram, INTENSITY);
}

int sinogram_attenuation(const field<double>& phantom, int projections,field<double>& sinogram){
    sinogram_base(phatom, projections, sinogram, ATTENUATION);
}