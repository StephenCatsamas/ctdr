
#include "field.h"
#include <cstdint>

#define PACK(c0, c1, c2, c3) \
    (((uint32_t)(uint8_t)(c0) << 24) | \
    ((uint32_t)(uint8_t)(c1) << 16) | \
    ((uint32_t)(uint8_t)(c2) << 8) | \
    ((uint32_t)(uint8_t)(c3)))

template<>
//maximum norm of two fields
double field<double>::norm(field<double>& a, field<double>& b){
    int h = a.height;
    int w = a.width;
    if(!(h == b.height && w == b.width)){exit(6);}

    double maxdiff = a[0][0];
    for(int i = 0; i< h; i++){
    for(int j = 0; j< w; j++){
        double diff = abs(a[i][j] - b[i][j]);
        maxdiff = std::max<double>(diff,maxdiff);     
    }
    }
    return maxdiff;
}

//reads in a png and gives a field of doubles 
//has two modes R_ONLY only takes only the R channel data
//RGBA reads in the RGBA data and interprets it as an
//IEEE754 float
field<double> png2doub(const char* fp, png_mode mode){
    auto [image,height,width] = decode_png(fp);

    auto fld = field<double>(height,width);

    for (int i = 0; i< height; i++){
    for (int j = 0; j< width; j++){
        int datindx = 4*(i*width + j);
        uint8_t r = image[datindx+0]; 
        uint8_t g = image[datindx+1]; 
        uint8_t b = image[datindx+2]; 
        uint8_t a = image[datindx+3]; 

        if (mode == R_ONLY){

        // fld[i][j] = ((int)r - 128); 
            fld[i][j] = r/255.0; 
        }else if (mode == RGBA){
            uint32_t col = PACK(r,g,b,a); 
            fld[i][j] = *(float*)&col;//interpret int data as float
        }else{
            out.log(ERR) << "invalid png reading mode\n";
        }
    }
    }
    return fld;
}

//reads in a png and gives a field of keepouts
field<int> png2int(const char* fp){
    auto [image,height,width] = decode_png(fp);

    auto bf = field<int>(height,width);

    for (int i = 0; i< height; i++){
    for (int j = 0; j< width; j++){
        int datindx = 4*(i*width + j);
        uint8_t r = image[datindx+0]; 
        uint8_t g = image[datindx+1]; 
        uint8_t b = image[datindx+2]; 
        uint8_t a = image[datindx+3]; 

        uint32_t col = PACK(r,g,b,a); 

        switch (col){
            case 0xFF0000FF: bf[i][j] = 1; break; //red is keepout
            default: bf[i][j] = 0;
        }
    }
    }
    return bf;
}

//takes in an int and writes it as RGBA data
void write_pixle(std::vector<unsigned char>& image, int p){
    uint8_t R = ((p >> 24) & 0xFF);
    uint8_t G = ((p >> 16) & 0xFF);
    uint8_t B = ((p >> 8) & 0xFF);
    uint8_t A = ((p >> 0) & 0xFF);

    image.push_back(R);
    image.push_back(G);
    image.push_back(B);
    image.push_back(A);
}



//saves double field to png
void doub2png(const char* fp, const field<double>& f, png_out_mode mode){
    auto image = std::vector<unsigned char>();
    int h = f.height;
    int w = f.width;

    double mn = f.min();
    double mx = f.max();

    uint8_t pixv;
    double vc;
    for(auto v : f.data){
        switch (mode){
            case RAW:
                vc = std::clamp(v, 0.0, 1.0);
                pixv = vc*255;
                write_pixle(image,PACK(pixv,pixv,pixv,0xFF));
                break;
            case SCALE:
                pixv = (v - mn)/(mx - mn)*255;
                write_pixle(image,PACK(pixv,255-pixv,pixv,0xFF));
                break;
            case SCALE_ZERO:
                pixv = (v/(std::max(abs(mx),abs(mn)))+0.5)*255;
                write_pixle(image,PACK(pixv,255-pixv,pixv,0xFF));
                break;
        }
        
    }

    encode_png(fp, image, w,h);
}

//unit testing
int test_field_interpolate(){
    auto A = field<double>(10,4);
    for(int i = 0; i < A.height; i++){
    for(int j = 0; j < A.width; j++){
        A[i][j] = i+100*j;
    }
    }

    std::cout << A << std::endl;

    auto B = A.subsample(1);

    std::cout << B << std::endl; 

    A.interpolate(B);

    std::cout << A << std::endl;

    return 1;
}