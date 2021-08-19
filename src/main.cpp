/**
 * Nikolaus Mayer, 2017 (mayern@cs.uni-freiburg.de)
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>

#include "CImg.h"
using namespace cimg_library;


int main(int argc, char** argv)
{
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <input> <output> [coeffs...]\n";
    return EXIT_FAILURE;
  }

  CImg<unsigned char> src;
  src.load(argv[1]);
  const int W{src.width()};
  const int H{src.height()};
  CImg<unsigned char> dst(W,H,1,src.spectrum());

  CImg<float> tmp(W,H,1,src.spectrum());
  CImg<unsigned char> tmp_display(W,H,1,src.spectrum());
  CImg<float> tmp_weights(W,H,1,1);
  CImg<float> map(W,H,1,2);
  CImg<float> inverse_map(W,H,1,2);

  const float cx{W/2.f};
  const float cy{H/2.f};

  const float c2{(argc >= 4) ? atof(argv[3]) : 2e-5};
  const float c4{(argc >= 5) ? atof(argv[4]) : 5e-7};
  const float c6{(argc >= 6) ? atof(argv[5]) : 1e-8};

  /// radial part of a radial-tangential camera distortion model
  auto f = [c2,c4,c6](float v) -> float {
    return v + c2*v*v + c4*v*v*v + c6*v*v*v*v;
  };

  /// approximately invert the distortion function
  const float threshold = 0.01;
  unsigned long int stepcount = 0;
  auto fi = [W,threshold,&f,&stepcount](float v) -> float {
    ++stepcount;
    int this_count = 0;
    float lower = 0.f;
    float upper = W;
    float guess = (lower+upper)/2;
    float fguess = f(guess);
    while (std::abs(fguess - v) > threshold and this_count < 100) {
      ++stepcount;
      ++this_count;
      if (fguess <= v)
        lower = guess;
      else
        upper = guess;
      guess = (lower+upper)/2;
      fguess = f(guess);
    }
    return guess;
  };

  for (int _y = 0; _y < H; ++_y) {
    for (int _x = 0; _x < W; ++_x) {
      float x = _x + .5f - cx;
      float y = _y + .5f - cy;
      float radius = std::sqrt(x*x + y*y);
      float nradius = f(radius);
      float nx = x*nradius/radius;
      float ny = y*nradius/radius;
      float _nx = nx - .5f + cx;
      float _ny = ny - .5f + cy;
      map(_x,_y,0,0) = _nx-_x;
      map(_x,_y,0,1) = _ny-_y;

      float iradius = fi(radius);
      nx = x*iradius/radius;
      ny = y*iradius/radius;
      _nx = nx - .5f + cx;
      _ny = ny - .5f + cy;
      inverse_map(_x,_y,0,0) = _nx-_x;
      inverse_map(_x,_y,0,1) = _ny-_y;
    }
  }

  std::cout << "Map inversion took "
            << static_cast<float>(stepcount)/(W*H)
            << " steps per pixel on average" << std::endl;

  /// Use inverse map to blit src onto tmp
  for (int _y = 0; _y < H; ++_y) {
    for (int _x = 0; _x < W; ++_x) {
      float x = _x+inverse_map(_x,_y,0,0);
      float y = _y+inverse_map(_x,_y,0,1);
      int xbase = static_cast<int>(std::trunc(x));
      int ybase = static_cast<int>(std::trunc(y));
      float right = x-xbase;
      float lower = y-ybase;
      if (xbase < 0 or xbase >= W-1 or ybase < 0 or ybase >= H-1)
        continue;
      tmp_weights(xbase  ,ybase  ) += (1.f-right)*(1.-lower);
      for (int _c = 0; _c < src.spectrum(); ++_c)
        tmp(xbase  ,ybase  ,0,_c) += src.linear_atXY(_x,_y,0,_c,0)*(1.f-right)*(1.-lower);

      tmp_weights(xbase+1,ybase  ) += (    right)*(1.-lower);
      for (int _c = 0; _c < src.spectrum(); ++_c)
        tmp(xbase+1,ybase  ,0,_c) += src.linear_atXY(_x,_y,0,_c,0)*(    right)*(1.-lower);

      tmp_weights(xbase  ,ybase+1) += (1.f-right)*(   lower);
      for (int _c = 0; _c < src.spectrum(); ++_c)
        tmp(xbase  ,ybase+1,0,_c) += src.linear_atXY(_x,_y,0,_c,0)*(1.f-right)*(   lower);

      tmp_weights(xbase+1,ybase+1) += (    right)*(   lower);
      for (int _c = 0; _c < src.spectrum(); ++_c)
        tmp(xbase+1,ybase+1,0,_c) += src.linear_atXY(_x,_y,0,_c,0)*(    right)*(   lower);
    }
  }

  /// Weight-normalize tmp
  for (int _y = 0; _y < H; ++_y) {
    for (int _x = 0; _x < W; ++_x) {
      for (int _c = 0; _c < src.spectrum(); ++_c) {
        tmp(_x,_y,0,_c) /= tmp_weights(_x,_y);
        tmp(_x,_y,0,_c) = std::max(0.f,std::min(tmp(_x,_y,0,_c),255.f));
        tmp_display(_x,_y,0,_c) = static_cast<unsigned char>(tmp(_x,_y,0,_c));
      }
    }
  }

  for (int _y = 0; _y < H; ++_y) {
    for (int _x = 0; _x < W; ++_x) {
      float x = _x+inverse_map(_x,_y,0,0);
      float y = _y+inverse_map(_x,_y,0,1);
      for (int _c = 0; _c < src.spectrum(); ++_c)
        dst(_x,_y,0,_c) = tmp_display.linear_atXY(x,y,0,_c,0);
    }
  }

  dst.save(argv[2]);
  tmp_display.save("tmp.ppm");

  /// Bye!
  return EXIT_SUCCESS;
}

