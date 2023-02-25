
#include "util.h"

#include "lodepng.h"

#include <ostream>
#include<tuple>
#include<vector>
#include <numbers>

double pi = std::numbers::pi_v<double>;
logger out = logger(INF, std::cout);

logger::logger(print_mode, std::ostream& output){
  printer = &output;
}

std::ostream& logger::log(print_mode mode){
  (*printer).clear(mode >= level
  ? std::ios_base::goodbit
  : std::ios_base::badbit);
  return *printer;
}

std::tuple<std::vector<unsigned char>,unsigned,unsigned> decode_png(const char* filename) {
  std::vector<unsigned char> image; //the raw pixels
  unsigned width, height;

  //decode
  unsigned error = lodepng::decode(image, width, height, filename);

  //if there's an error, display it
  if(error) {out.log(ERR) << "decoder error " << error << " on file " << filename << ": " << lodepng_error_text(error) << std::endl; exit(3);}

  //the pixels are now in the vector "image", 4 bytes per pixel, ordered RGBARGBA..., use it as texture, draw it, ...
  return {image,height,width};
}

void encode_png(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
  unsigned error = lodepng::encode(filename, image, width, height);

  if(error) out.log(ERR) << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}