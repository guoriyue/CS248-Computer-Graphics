#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CS248 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Advanced Task
  // Implement mipmap for trilinear filtering

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 4: Implement nearest neighbour interpolation
  
  if (level < 0 || level > tex.mipmap.size()){
    // return magenta for invalid level
    return Color(1,0,1,1);
  }

  int width = tex.mipmap[level].width;
  int height = tex.mipmap[level].height;

  int x = (int) round(u * width - 0.5);
  int y = (int) round(v * height - 0.5);

  Color texuture_color;
  texuture_color.r = tex.mipmap[level].texels[4 * (x + y * width)] / 255.0;
  texuture_color.g = tex.mipmap[level].texels[4 * (x + y * width) + 1] / 255.0;
  texuture_color.b = tex.mipmap[level].texels[4 * (x + y * width) + 2] / 255.0;
  texuture_color.a = tex.mipmap[level].texels[4 * (x + y * width) + 3] / 255.0;

  return texuture_color;

}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 4: Implement bilinear filtering

  if (level < 0 || level > tex.mipmap.size()){
    // return magenta for invalid level
    return Color(1,0,1,1);
  }
  int width = tex.mipmap[level].width;
  int height = tex.mipmap[level].height;

  float u_texture = u * width - 0.5;
  float v_texture = v * height - 0.5;

  if (u_texture < 0) u_texture = 0;
  if (v_texture < 0) v_texture = 0;
  if (u_texture > width-1) u_texture = width-1;
  if (v_texture > height-1) v_texture = height-1;

  int u0 = (int)u_texture;
  int u1 = u0 + 1;
  int v0 = (int)v_texture;
  int v1 = v0 + 1;

  Color color_00;
  color_00.r = tex.mipmap[level].texels[4 * (u0 + v0 * width)] / 255.0;
  color_00.g = tex.mipmap[level].texels[4 * (u0 + v0 * width) + 1] / 255.0;
  color_00.b = tex.mipmap[level].texels[4 * (u0 + v0 * width) + 2] / 255.0;
  color_00.a = tex.mipmap[level].texels[4 * (u0 + v0 * width) + 3] / 255.0;

  Color color_10;
  color_10.r = tex.mipmap[level].texels[4 * (u1 + v0 * width)] / 255.0;
  color_10.g = tex.mipmap[level].texels[4 * (u1 + v0 * width) + 1] / 255.0;
  color_10.b = tex.mipmap[level].texels[4 * (u1 + v0 * width) + 2] / 255.0;
  color_10.a = tex.mipmap[level].texels[4 * (u1 + v0 * width) + 3] / 255.0;

  Color color_01;
  color_01.r = tex.mipmap[level].texels[4 * (u0 + v1 * width)] / 255.0;
  color_01.g = tex.mipmap[level].texels[4 * (u0 + v1 * width) + 1] / 255.0;
  color_01.b = tex.mipmap[level].texels[4 * (u0 + v1 * width) + 2] / 255.0;
  color_01.a = tex.mipmap[level].texels[4 * (u0 + v1 * width) + 3] / 255.0;

  Color color_11;
  color_11.r = tex.mipmap[level].texels[4 * (u1 + v1 * width)] / 255.0;
  color_11.g = tex.mipmap[level].texels[4 * (u1 + v1 * width) + 1] / 255.0;
  color_11.b = tex.mipmap[level].texels[4 * (u1 + v1 * width) + 2] / 255.0;
  color_11.a = tex.mipmap[level].texels[4 * (u1 + v1 * width) + 3] / 255.0;

  Color color_u0;
  color_u0.r = color_00.r + ((u_texture - u0)/ (u1 - u0)) * (color_10.r - color_00.r);
  color_u0.g = color_00.g + ((u_texture - u0)/ (u1 - u0)) * (color_10.g - color_00.g);
  color_u0.b = color_00.b + ((u_texture - u0)/ (u1 - u0)) * (color_10.b - color_00.b);
  color_u0.a = color_00.a + ((u_texture - u0)/ (u1 - u0)) * (color_10.a - color_00.a);

  Color color_u1;
  color_u1.r = color_01.r + ((u_texture - u0)/ (u1 - u0)) * (color_11.r - color_01.r);
  color_u1.g = color_01.g + ((u_texture - u0)/ (u1 - u0)) * (color_11.g - color_01.g);
  color_u1.b = color_01.b + ((u_texture - u0)/ (u1 - u0)) * (color_11.b - color_01.b);
  color_u1.a = color_01.a + ((u_texture - u0)/ (u1 - u0)) * (color_11.a - color_01.a);

  Color color_texure;
  color_texure.r = color_u0.r + ((v_texture - v0)/ (v1 - v0)) * (color_u1.r - color_u0.r);
  color_texure.g = color_u0.g + ((v_texture - v0)/ (v1 - v0)) * (color_u1.g - color_u0.g);
  color_texure.b = color_u0.b + ((v_texture - v0)/ (v1 - v0)) * (color_u1.b - color_u0.b);
  color_texure.a = color_u0.a + ((v_texture - v0)/ (v1 - v0)) * (color_u1.a - color_u0.a);

  return color_texure;

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Advanced Task
  // Implement trilinear filtering

  // return magenta for invalid level
  return Color(1,0,1,1);

}

} // namespace CS248
