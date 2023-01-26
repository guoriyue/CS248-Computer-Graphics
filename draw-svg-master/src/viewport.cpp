#include "viewport.h"

#include "CS248.h"

namespace CS248 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 3 (part 2): 
  // Set svg to normalized device coordinate transformation. 
  // Your input arguments are defined as SVG canvans coordinates.

  this->x = x;
  this->y = y;
  this->span = span; 
  
  Matrix3x3 m1 = Matrix3x3::identity();
  m1(0,0) = 1.0/(2*span); m1(0,1) = 0; m1(0,2) = 0;
  m1(1,0) = 0; m1(1,1) = 1.0/(2*span); m1(1,2) = 0;
  m1(2,0) = 0; m1(2,1) = 0; m1(2,2) = 1;

  Matrix3x3 m2 = Matrix3x3::identity();
  m2(0,0) = 1; m2(0,1) = 0; m2(0,2) = span - x;
  m2(1,0) = 0; m2(1,1) = 1; m2(1,2) = span - y;
  m2(2,0) = 0; m2(2,1) = 0; m2(2,2) = 1;
  
  
  set_canvas_to_norm(m1 * m2);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CS248
