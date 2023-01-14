#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CS248 {


// Implements SoftwareRenderer //

// fill a sample location with color
void SoftwareRendererImp::fill_sample(int sx, int sy, const Color &color) {
  // Task 2: implement this function
  // printf("fill_sample");
  if (sx < 0 || sx >= sample_width) return;
  if (sy < 0 || sy >= sample_height) return;

  Color pixel_color;
  float inv255 = 1.0 / 255.0;
  pixel_color.r = sample_buffer[4 * (sx + sy * sample_width)] * inv255;
  pixel_color.g = sample_buffer[4 * (sx + sy * sample_width) + 1] * inv255;
  pixel_color.b = sample_buffer[4 * (sx + sy * sample_width) + 2] * inv255;
  pixel_color.a = sample_buffer[4 * (sx + sy * sample_width) + 3] * inv255;

  pixel_color = ref->alpha_blending_helper(pixel_color, color);

  sample_buffer[4 * (sx + sy * sample_width)] = (uint8_t)(pixel_color.r * 255);
  sample_buffer[4 * (sx + sy * sample_width) + 1] = (uint8_t)(pixel_color.g * 255);
  sample_buffer[4 * (sx + sy * sample_width) + 2] = (uint8_t)(pixel_color.b * 255);
  sample_buffer[4 * (sx + sy * sample_width) + 3] = (uint8_t)(pixel_color.a * 255);
}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

	// // Task 2: Re-implement this function

	// // check bounds
	// if (x < 0 || x >= width) return;
	// if (y < 0 || y >= height) return;

	// Color pixel_color;
	// float inv255 = 1.0 / 255.0;
	// pixel_color.r = pixel_buffer[4 * (x + y * width)] * inv255;
	// pixel_color.g = pixel_buffer[4 * (x + y * width) + 1] * inv255;
	// pixel_color.b = pixel_buffer[4 * (x + y * width) + 2] * inv255;
	// pixel_color.a = pixel_buffer[4 * (x + y * width) + 3] * inv255;

	// pixel_color = ref->alpha_blending_helper(pixel_color, color);

	// pixel_buffer[4 * (x + y * width)] = (uint8_t)(pixel_color.r * 255);
	// pixel_buffer[4 * (x + y * width) + 1] = (uint8_t)(pixel_color.g * 255);
	// pixel_buffer[4 * (x + y * width) + 2] = (uint8_t)(pixel_color.b * 255);
	// pixel_buffer[4 * (x + y * width) + 3] = (uint8_t)(pixel_color.a * 255);

  int sx = sample_rate*x;
  int sy = sample_rate*y;
  if (sx < 0 || sx >= sample_width) return;
	if (sy < 0 || sy >= sample_height) return;

  for(int i=sx; i<sx+sample_rate; i++){
    for(int j=sy; j<sy+sample_rate; j++){
      fill_sample(i, j, color);
    }
  }
}

void SoftwareRendererImp::draw_svg( SVG& svg ) {
  printf("draw_svg");
  this->sample_buffer = new unsigned char[this->sample_width*this->sample_height];
  // set top level transformation
  transformation = canvas_to_screen;

  // canvas outline
  Vector2D a = transform(Vector2D(0, 0)); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width, 0)); b.x++; b.y--;
  Vector2D c = transform(Vector2D(0, svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width, svg.height)); d.x++; d.y++;

  svg_bbox_top_left = Vector2D(a.x+1, a.y+1);
  svg_bbox_bottom_right = Vector2D(d.x-1, d.y-1);

  // draw all elements
  for (size_t i = 0; i < svg.elements.size(); ++i) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to pixel buffer
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {
  printf("set_sample_rate");
  // Task 2: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;

  this->sample_width = width*sample_rate;
  this->sample_height = height*sample_rate;
  // this->sample_buffer = (unsigned char*) malloc(sizeof(unsigned char));
  this->sample_buffer = new unsigned char[this->sample_width*this->sample_height];
}

void SoftwareRendererImp::set_pixel_buffer( unsigned char* pixel_buffer,
                                             size_t width, size_t height ) {
                                              
  printf("set_pixel_buffer");
  // Task 2: 
  // You may want to modify this for supersampling support
  this->pixel_buffer = pixel_buffer;
  this->width = width;
  this->height = height;

  this->sample_width = width*sample_rate;
  this->sample_height = height*sample_rate;
  // this->sample_buffer = (unsigned char*) malloc(sizeof(unsigned char));
  this->sample_buffer = new unsigned char[this->sample_width*this->sample_height];
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

	// Task 3 (part 1):
	// Modify this to implement the transformation stack

	switch (element->type) {
	case POINT:
		draw_point(static_cast<Point&>(*element));
		break;
	case LINE:
		draw_line(static_cast<Line&>(*element));
		break;
	case POLYLINE:
		draw_polyline(static_cast<Polyline&>(*element));
		break;
	case RECT:
		draw_rect(static_cast<Rect&>(*element));
		break;
	case POLYGON:
		draw_polygon(static_cast<Polygon&>(*element));
		break;
	case ELLIPSE:
		draw_ellipse(static_cast<Ellipse&>(*element));
		break;
	case IMAGE:
		draw_image(static_cast<Image&>(*element));
		break;
	case GROUP:
		draw_group(static_cast<Group&>(*element));
		break;
	default:
		break;
	}

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Advanced Task
  // Implement ellipse rasterization

}

void SoftwareRendererImp::draw_image( Image& image ) {

  // Advanced Task
  // Render image element with rotation

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;

  // fill sample - NOT doing alpha blending!
  // TODO: Call fill_pixel here to run alpha blending
  pixel_buffer[4 * (sx + sy * width)] = (uint8_t)(color.r * 255);
  pixel_buffer[4 * (sx + sy * width) + 1] = (uint8_t)(color.g * 255);
  pixel_buffer[4 * (sx + sy * width) + 2] = (uint8_t)(color.b * 255);
  pixel_buffer[4 * (sx + sy * width) + 3] = (uint8_t)(color.a * 255);

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 0: 
  // Implement Bresenham's algorithm (delete the line below and implement your own)
  ref->rasterize_line_helper(x0, y0, x1, y1, width, height, color, this);

  // Advanced Task
  // Drawing Smooth Lines with Line Width
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  printf("rasterize_triangle");
  // Task 1: 
  // Implement triangle rasterization
  
  int min_x = floor(min(min(x0, x1), x2));
  int max_x = ceil(max(max(x0, x1), x2));
  int min_y = floor(min(min(y0, y1), y2));
  int max_y = ceil(max(max(y0, y1), y2));

  if (min_x < 0) min_x=0;
  if (min_y < 0) min_y=0;
  if (max_x > width) max_x = ceil(width);
  if (max_y > height) max_y = ceil(height);
  
  float b[3] = {x1-x0, x2-x1, x0-x2};
  float a[3] = {y1-y0, y2-y1, y0-y2};
  float c[3] = {y0*(x1-x0)-x0*(y1-y0), y1*(x2-x1)-x1*(y2-y1), y2*(x0-x2)-x2*(y0-y2)};
  bool anti_clockwise = x0*y1+y0*x2+x1*y2-x2*y1-y2*x0-x1*y0>0 ? true:false;
  // https://math.stackexchange.com/questions/1324179/how-to-tell-if-3-connected-points-are-connected-clockwise-or-counter-clockwise
  // larger than 0, anti clockwise order

  for(int i=min_x; i<=max_x; i++){
    for(int j=min_y; j<=max_y; j++){
      float sample_x=i+0.5;
      float sample_y=j+0.5;
      float l[3]={a[0]*sample_x-b[0]*sample_y+c[0], a[1]*sample_x-b[1]*sample_y+c[1], a[2]*sample_x-b[2]*sample_y+c[2]};
      if((anti_clockwise && l[0]<=0 && l[1]<=0 && l[2]<=0) || (!anti_clockwise && l[0]>=0 && l[1]>=0 && l[2]>=0)){
        fill_pixel(sample_x, sample_y, color);
      }
    }
  }
  // 1 pixel different for stanford
  // 0 for traingles
  // 1 for cube
  // 0 for kitten


  // Advanced Task
  // Implementing Triangle Edge Rules

}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 4: 
  // Implement image rasterization

}

// resolve samples to pixel buffer
void SoftwareRendererImp::resolve( void ) {

  // Task 2: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 2".
  printf("resolve");
  for(int i=0; i<width; i++){
    for(int j=0; j<height; j++){
      float r=0.0; float g=0.0; float b=0.0; float a=0.0;
      for(int si=0; si<sample_rate; si++){
        for(int sj=0; sj<sample_rate; sj++){
          int sx=sample_rate*i+si;
          int sy=sample_rate*j+sj;
          r+=sample_buffer[4*(sx+sy*sample_width)];
          g+=sample_buffer[4*(sx+sy*sample_width)+1];
          b+=sample_buffer[4*(sx+sy*sample_width)+2];
          a+=sample_buffer[4*(sx+sy*sample_width)+3];
        }
      }
      r/=(sample_rate*sample_rate);
      g/=(sample_rate*sample_rate);
      b/=(sample_rate*sample_rate);
      a/=(sample_rate*sample_rate);
      pixel_buffer[4*(i+j*width)]=r;
      pixel_buffer[4*(i+j*width)+1]=g;
      pixel_buffer[4*(i+j*width)+2]=b;
      pixel_buffer[4*(i+j*width)+3]=a;
    }
  }
  // free(sample_buffer);
  return;
}

Color SoftwareRendererImp::alpha_blending(Color pixel_color, Color color)
{
  // Task 5
  // Implement alpha compositing
  return pixel_color;
}


} // namespace CS248
