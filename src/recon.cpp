
#include "recon.h"
#include "util.h"
#include "field.h"
#include "pocketfft.h"

using namespace pocketfft;

int sinogram2projection(const field<double>& sinogram, const double angle, std::vector<double>& projection){
    
    projection.resize(sinogram.width);  

    int angle_index = (int)std::round(angle/(2.0*pi)*sinogram.height) %sinogram.height;//if we always round down then we will get some rotation errors
    auto sino_iter = sinogram.data.begin();

    std::copy(  sino_iter+(sinogram.width* angle_index   ), 
                sino_iter+(sinogram.width*(angle_index+1)), 
                projection.begin());
    
    return 1;
}

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

int back_project(const std::vector<double>& projection, double angle, field<double>& back_projection){
    
    for (int j = 0; j< back_projection.width; j++){
    for (int i = 0; i< back_projection.height; i++){
        back_projection[i][j] = projection[j];
    }
    } 
    back_projection.rot(-angle);
    
    return 1;
}




int recon_bp(const field<double>& sinogram, const int n_proj, field<double>& tomogram){
    
    const double angle_step = 2.0*pi/n_proj;
    
    tomogram = field<double>(sinogram.width, sinogram.width, 0.0);
    auto projection = std::vector<double>(sinogram.width);
    auto back_projection = field<double>(sinogram.width,sinogram.width);
    for(int i = 0; i < n_proj; i++){
        double angle = i * angle_step;
        
        sinogram2projection(sinogram,angle,projection);
        
        projection *= 8.0*(0.5/(n_proj*back_projection.height));
        back_project(projection, angle, back_projection);

        tomogram += back_projection ;
       
    }
    return 1; 
}

int recon_fbp(const field<double>& sinogram, const int n_proj, field<double>& tomogram){
    
    const double angle_step = 2.0*pi/n_proj;   
    
    tomogram = field<double>(sinogram.width, sinogram.width, 0.0);
    std::vector<double> projection;
    auto back_projection = field<double>(tomogram.height,tomogram.width);
    for(double angle = 0.0; angle < 2.0*pi; angle += angle_step){
        // out.log(DBG) << angle << std::endl;
        sinogram2projection(sinogram,angle,projection);
        
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
        for(int i = 0; i < f_projection.size(); i++){
            double weight = 2*pi*i;
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
        projection *= (0.5/(n_proj*back_projection.height));        
        back_project(projection, angle, back_projection);
        tomogram += back_projection;
       
    }
    return 1; 
}


int recon_dfi(const field<double>& sinogram, const int n_proj, const int pf, field<double>& tomogram){
    
    const double angle_step = 2.0*pi/n_proj;   
    auto f_polar_proj = field<std::complex<double>>(n_proj,pf*tomogram.width/2 + 1);
    auto f_tomogram = field<std::complex<double>>(tomogram.height,tomogram.width, 0.0);
    std::vector<double> projection;
    tomogram = field<double>(sinogram.width, sinogram.width, 0.0);
    //polar fouier space
    for(int i = 0; i < n_proj; i++){
        double angle = i * angle_step;
        // out.log(INF) << angle << std::endl;
        
        sinogram2projection(sinogram,angle,projection);
        
        projection.resize(pf*projection.size(), 0.0);
        std::rotate(projection.begin(),projection.begin()+sinogram.width/2,projection.end());//rotation very importaint to make phase changes slowly
        
        shape_t shape = {projection.size()};
        stride_t stride_r = {sizeof(double)};
        stride_t stride_c = {sizeof(std::complex<double>)};
        shape_t axes = {0};

        auto f_projection = std::vector<std::complex<double>>(shape[0]/2 + 1);
        double scale_factor = 0.5/(sinogram.width*sinogram.width);

        r2c(shape,
            stride_r, 
            stride_c, 
            axes,
            FORWARD, 
            projection.data(), 
            f_projection.data(), 
            scale_factor);
        
        std::copy(f_projection.begin(), f_projection.end(), 
                  f_polar_proj.data.begin()+(f_polar_proj.width*i));       
    }

    pol2cart(f_polar_proj, f_tomogram);

    
    shape_t shape_2d = {(size_t)f_tomogram.width, 
                        (size_t)f_tomogram.height};
    stride_t stride_2d = {sizeof(std::complex<double>),
                          sizeof(std::complex<double>)*f_tomogram.width};
    shape_t axes_2d = {0,1};   
               
    c2c(shape_2d,
        stride_2d, 
        stride_2d, 
        axes_2d,
        BACKWARD, 
        f_tomogram.data.data(), 
        f_tomogram.data.data(), 
        1.0);
    
    for(int i = 0; i< f_tomogram.height; i++){
    for(int j = 0; j< f_tomogram.width; j++){
        int ic = (i+f_tomogram.height/2)%f_tomogram.height;
        int jc = (j+f_tomogram.width/2)%f_tomogram.width;
        tomogram[i][j] = 2.0*std::abs(f_tomogram[ic][jc]);
    }
    }
    {
    // out.log(INF) << "saving" << std::endl;
    // field<double> real_field;
    // real_field = abs(f_polar_proj);
    // doub2png("data/sinogram_abs.png", real_field, SCALE);
    // real_field = arg(f_polar_proj);
    // doub2png("data/sinogram_arg.png", real_field, SCALE_ZERO);
    
    // real_field = abs(f_tomogram);
    // doub2png("data/f_tomogram_abs.png", real_field, SCALE);
    // real_field = arg(f_tomogram);
    // doub2png("data/f_tomogram_arg.png", real_field, SCALE_ZERO);
    }
     
    return 1;
}

int recon_art(const field<double>& sinogram, const int n_proj, field<double>& tomogram){
    
    tomogram = field<double>(sinogram.width, sinogram.width, 0.0);
    const double angle_step = 2.0*pi/n_proj; 
    auto tomo_delta = field<double>(tomogram.height,tomogram.width);
    std::vector<double> p;
    std::vector<double> q;

    for(int i = 0; i < n_proj; i++){
        size_t memememe = sizeof(i);
        size_t fhfhfh = sizeof(n_proj);
        double angle = (1500450271LL * i)%n_proj * angle_step;//use a big prime to do the correction in a semi random order 
            
        sinogram2projection(sinogram,angle,p);
        project(tomogram, angle, q);
        
        auto correction = p-q;
        correction *= (0.8/(tomo_delta.height));
        back_project(correction, angle, tomo_delta);
        tomogram += tomo_delta;
            
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

