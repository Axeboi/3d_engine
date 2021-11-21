#pragma once

#include <iostream>
#include <algorithm>
//#include <cmath>
#include <cassert>
#include <cstdlib>

#include "render_data_structures.h"

const double PI = 3.14159265358979323846;

// Forward declarations
void draw_triangle(std::vector<std::vector<ScreenBuffer>> &, size_t, size_t, Vec4&, Vec4& , Vec4& , float );
bool float_equal(float, float);

void load_scene(TriangleBuffer &buffer)
{
    // TEST TRIANGLE
   /* triangle_buffer_add_triangle(buffer, {
        -0.5, -0.5, 4.0, 1.0,
         0.0,  0.5, 4.0, 1.0,
         0.5, -0.5, 4.0, 1.0,
    });*/
  triangle_buffer_add_triangle(buffer, {
    -0.5, -0.5, -0.5, 1.0,
    -0.5, 0.5, -0.5, 1.0,
     0.5, 0.5, -0.5, 1.0
  });
  
  triangle_buffer_add_triangle(buffer, {
    -0.5, -0.5, -0.5, 1.0,
    0.5, 0.5, -0.5, 1.0,
    0.5, -0.5, -0.5, 1.0,
  });

  triangle_buffer_add_triangle(buffer, {
    0.5, -0.5, -0.5, 1.0,
    0.5, 0.5, -0.5, 1.0,
    0.5, 0.5, 0.5, 1.0
  });

  triangle_buffer_add_triangle(buffer, {
    0.5, -0.5, -0.5, 1.0,
	  0.5, 0.5, 0.5, 1.0,
	  0.5, -0.5, 0.5, 1.0
  });

  triangle_buffer_add_triangle(buffer, {
    0.5, -0.5, 0.5, 1.0,
	  0.5, 0.5, 0.5, 1.0,
	  -0.5, 0.5, 0.5, 1.0
  });

  triangle_buffer_add_triangle(buffer, {
    0.5, -0.5, 0.5, 1.0,
	  -0.5, 0.5, 0.5, 1.0,
	  -0.5, -0.5, 0.5, 1.0
  });
  
  triangle_buffer_add_triangle(buffer, {
    -0.5, -0.5, 0.5, 1.0,
	  -0.5, 0.5, 0.5, 1.0,
	  -0.5, 0.5, -0.5, 1.0
  });
  
  triangle_buffer_add_triangle(buffer, {
    -0.5, -0.5, 0.5, 1.0,
	  -0.5, 0.5, -0.5, 1.0,
	  -0.5, -0.5, -0.5, 1.0
  });
  
  triangle_buffer_add_triangle(buffer, {
    -0.5, 0.5, -0.5, 1.0,
	  -0.5, 0.5, 0.5, 1.0,
	  0.5, 0.5, 0.5, 1.0
  });
  
  
  triangle_buffer_add_triangle(buffer, {
    -0.5, 0.5, -0.5, 1.0,
	  0.5, 0.5, 0.5, 1.0,
	  0.5, 0.5, -0.5, 1.0
  });
 
  triangle_buffer_add_triangle(buffer, {
   	0.5, -0.5, 0.5, 1.0,
	  -0.5, -0.5, 0.5, 1.0,
	  -0.5, -0.5, -0.5, 1.0
  });
	
  triangle_buffer_add_triangle(buffer, {
    0.5, -0.5, 0.5, 1.0,
	  -0.5, -0.5, -0.5, 1.0,
	  0.5, -0.5, -0.5, 1.0
  });
};

void transform_vectors(TriangleBuffer &triangles, Vec4 &move)
{   
    float theta = 3.1415;
    Matrix4 x_rot {
      1, 0, 0, 0,
      0, std::cos(theta/2), -std::sin(theta/2), 0,
      0, std::sin(theta/2), std::cos(theta/2), 0,
      0, 0, 0, 1,
    };

    Matrix4 z_rot {
      std::cos(theta + theta/4), -std::sin(theta + theta/4), 0, 0, 
      std::sin(theta + theta/4), std::cos(theta + theta/4), 0, 0,
      0, 0, 1, 0, 
      0, 0, 0, 1,
    };

    Matrix4 y_rot {
      std::cos(theta/2), 0, std::sin(theta/2), 0,
      0, 1, 0, 0,
      -std::sin(theta/2), 0, std::cos(theta/2), 0,
      0, 0, 0, 1,
    };

    Matrix4 x_rot_inverse {
      1, 0, 0, 0,
      0, -std::cos(theta/2), std::sin(theta/2), 0,
      0, -std::sin(theta/2), -std::cos(theta/2), 0,
      0, 0, 0, 1,
    };

     Matrix4 z_rot_inverse {
      -std::cos(theta/2), std::sin(theta/2), 0, 0,
      -std::sin(theta/2), -std::cos(theta/2), 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
    };

    Matrix4 y_rot_inverse {
      -std::cos(theta/2), 0, -std::sin(theta/2), 0,
      0, 1, 0, 0,
      -std::sin(theta/2), 0, -std::cos(theta/2), 0,
      0, 0, 0, 1,
    };

    Matrix4 trans {
      1, 0, 0, 0,
      0, 1, 0, 0, 
      0, 0, 1, 0,
      0, 0, 400, 1
    };

    Matrix4 tmp;
    Matrix4 rotate_and_move;
    Vec4 tmp2;
    tmp = GeometryCalc::matrix_matrix_mul(z_rot, x_rot);
    rotate_and_move = tmp;
    tmp = GeometryCalc::matrix_matrix_mul(rotate_and_move, trans);
    rotate_and_move = tmp;
 

    tmp = GeometryCalc::matrix_matrix_mul(z_rot, x_rot);
    //tmp = GeometryCalc::vector_add(triangles.tris[i].tri[j], move);

    for (uint32_t i = 0; i < triangles.size; i++)
        for (uint32_t j = 0; j < 3; j++)
        {
            //tmp = GeometryCalc::vector_add(triangles.tris[i].tri[j], move);
            //tmp = GeometryCalc::matrix_matrix_mul(x_rot, trans);
            tmp2 = GeometryCalc::vec_matrix_mul(triangles.tris[i].tri[j], rotate_and_move);
            triangles.tris[i].tri[j] = tmp2;
            Vec4 n = GeometryCalc::vec_matrix_mul(triangles.tris[i].normal, tmp);
            triangles.tris[i].normal = n;
        }
};

void get_camera_transform(TriangleBuffer &buffer, Matrix4 &transform, Vec4 &camera_location, Vec4 &camera_direction)
{
  // ANGLES FOR CAMERA TRANSFORM AND ROTATION TO -Z AXIS:
    // First: calculate theta and rotate around the y axis
    // Second: calculate phi and rotate arounyd the x axis
    // Now the camera will point to -z axis

    // As long as the camera isn't rotating, a single vector should suffice to calculate camera
    // Should we want to add rotation, there needs to be an up and right vector also
    // example link for explanation: https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/lookat-function

    //Matrix4 transform;
    Matrix4 rotation;
    Matrix4 res;

    transform.m_data[3][0] += camera_location.x;
    transform.m_data[3][1] += camera_location.y;
    transform.m_data[3][2] += camera_location.z;
    
	/*double theta = -std::atan2(camera_direction.x, camera_direction.z);
	double phi = std::atan2(camera_direction.y, std::sqrt((camera_direction.x * camera_direction.x) + (camera_direction.y * camera_direction.y)));
    
    Matrix4 y_rotation {
          std::cos(theta),  0, -std::sin(theta), 0,
          0, 1.0, 0, 0,
          std::sin(theta), 0, std::cos(theta), 0,
          0, 0, 0, 1.0
    };

    Matrix4 x_rotation {
        1.0, 0, 0, 0,
        0, std::cos(phi),  std::sin(phi), 0,
        0, -std::sin(phi), std::cos(phi), 0,
        0, 0, 0, 1.0
    };

    GeometryCalc::matrix_row_to_column(x_rotation);
    rotation = GeometryCalc::matrix_matrix_mul(transform, y_rotation);

    GeometryCalc::matrix_row_to_column(rotation);
    res = GeometryCalc::matrix_matrix_mul(rotation, x_rotation);*/

    // for (uint32_t t = 0; t < buffer.size; t++)
    //     for (uint32_t v = 0; v < 3; v++)
    //     {
    //         Vec4 tmp = GeometryCalc::vec_matrix_mul(buffer.tris[t].tri[v], transform);
    //         buffer.tris[t].tri[v].x = tmp.x;
    //         buffer.tris[t].tri[v].y = tmp.y;
    //         buffer.tris[t].tri[v].z = tmp.z;
    //         buffer.tris[t].tri[v].w = tmp.w;
    //     }

  // translate camera to origin
  // rotate camera to -z axis
  // set transform
};

Matrix4 get_projection_transform(Matrix4 &transform, float f, float n, float fov)
{
  float magic_nr1 = (f + n) / (f - n);
  float magic_nr2 = (2 * f * n) / (f - n);

  Matrix4 projection {
    (std::cos(fov/2) / std::sin(fov/2)), 0,                                   0,          0,
    0,                                   (std::cos(fov/2) / std::sin(fov/2)), 0,          0,
    0,                                   0,                                   magic_nr1,  magic_nr2,
    0,                                   0,                                   -1,  0,
  };

  Matrix4 tmp = GeometryCalc::matrix_matrix_mul(transform, projection);
  transform = tmp;
  return transform;
  //return projection;
};

void apply_world_to_projection_transform(TriangleBuffer& triangles, Matrix4 transform)
{
    Vec4 tmp;
    for (size_t i = 0; i < triangles.size; i++)
        for (size_t v = 0; v < 3; v++)
        {
            tmp = GeometryCalc::vec_matrix_mul(triangles.tris[i].tri[v], transform);
            triangles.tris[i].tri[v] = tmp;
        }
}
     

void apply_light(TriangleBuffer &triangles, Vec4 light_source_normal) // 
{
  GeometryCalc::normalize(light_source_normal);
  Matrix4 camera_matrix {
    -1.0, 0, 0, 0,
    0, 1.0, 0, 0,
    0, 0, -1.0, 0,
    light_source_normal.x, light_source_normal.y, light_source_normal.z, 1,
  };
  
  for (size_t i = 0; i < triangles.size; i++)
  {
    Vec4 light_direction { 0.0, 0.0, -1.0, 1.0 };
    Vec4 point_location = GeometryCalc::vec_matrix_mul(triangles.tris[i].tri[0], camera_matrix);
    // if(GeometryCalc::vector_dot(point_location, light_direction) > 0)   // If triangle is visible to the camera, this should be a positive number
    // {

      Vec4 normal = triangles.tris[i].normal;
      GeometryCalc::normalize(normal);
      GeometryCalc::normalize(point_location);

      triangles.tris[i].shadow_interpolated = std::max(GeometryCalc::vector_cos(normal, point_location), 0.0f);
    // }
  }
};

void rasterize(const TriangleBuffer &triangles, std::vector<std::vector<ScreenBuffer>> &screen_buffer, size_t width, size_t height)
{
  for (size_t i = 0; i < triangles.size; i++)
  {
    // Triangles with w value less then 0 is behind the camera in the 4th dimension and will appear in front of the camera but mirrored
    // Here we need to do some clipping to make sure that triangle points close to the camera are clamped to w = 0;
    if(triangles.tris[i].tri[0].w > 0 && triangles.tris[i].tri[1].w > 0 && triangles.tris[i].tri[2].w > 0)
      continue;

    for (size_t p = 0; p < 3; p++)
    {
      triangles.tris[i].tri[p].x /= triangles.tris[i].tri[p].w;
      triangles.tris[i].tri[p].y /= triangles.tris[i].tri[p].w;
      triangles.tris[i].tri[p].z /= triangles.tris[i].tri[p].w;
    }

    draw_triangle(screen_buffer, width, height, triangles.tris[i].tri[0], triangles.tris[i].tri[1], triangles.tris[i].tri[2], triangles.tris[i].shadow_interpolated);
  }
};

void setImagePixel(std::vector<std::vector<ScreenBuffer>> &screen_buffer, const Point &p, int w0, int w1, int w2, float z0, float z1, float z2, float col)
{
  float fw0 = (float) w0 / (float) (w0 + w1 + w2);
  float fw1 = (float) w1 / (float) ((w0 + w1 + w2));
  float fw2 = (float) w2 / (float) ((w0 + w1 + w2));

  float z = z0 + (fw1 * (z1-z0)) + (fw2 * (z2 - z0));

  if (screen_buffer[p.y][p.x].z < z)
  {
    screen_buffer[p.y][p.x].z = z;
    // screen_buffer[p.y][p.x].r = fw0;
    // screen_buffer[p.y][p.x].g = fw1;
    // screen_buffer[p.y][p.x].b = fw2;
     screen_buffer[p.y][p.x].r = (unsigned int) (col * 255.0f);
     screen_buffer[p.y][p.x].g = (unsigned int) (col * 255.0f);
     screen_buffer[p.y][p.x].b = (unsigned int) (col * 255.0f);
    //do moar stuff
  }
}

bool is_top_left(const Point& a, const Point& b)
{
  
   return !(a.x == b.x || a.y < b.y);
}

// Calculate barycentric coordinates, one line at a time
// Note: swapping b and c will change the winding of triangle
int orient2d(const Point& a, const Point& b, const Point& c)
{
  return ((b.x - a.x) * (c.y - a.y)) - ((b.y - a.y) * (c.x - a.x));
}

int min3(int a, int b, int c) { return std::min(std::min(a, b), c); };
int max3(int a, int b, int c) { return std::max(std::max(a, b), c); };

void draw_triangle(std::vector<std::vector<ScreenBuffer>> &screen_buffer, size_t screen_width, size_t screen_height, Vec4& input_p0, Vec4& input_p1, Vec4& input_p2, float col)
{
  Point p0, p1, p2;

  //convert from image space (-1 to 1 x, -1 to 1 y) to screen space
  p0.x = (input_p0.x + 1) * (0.5 * (screen_width -1) );
  p0.y = (input_p0.y + 1) * (0.5 * (screen_height -1));
  p0.z = input_p0.z;

  p1.x = (input_p1.x + 1) * (0.5 * (screen_width -1));
  p1.y = (input_p1.y + 1) * (0.5 * (screen_height -1));
  p1.z = input_p1.z;

  p2.x = (input_p2.x + 1) * (0.5 * (screen_width -1));
  p2.y = (input_p2.y + 1) * (0.5 * (screen_height -1));
  p2.z = input_p2.z;
  
  // AABB for triangle
  int minX = min3(p0.x, p1.x, p2.x);
  int minY = min3(p0.y, p1.y, p2.y);
  int maxX = max3(p0.x, p1.x, p2.x);
  int maxY = max3(p0.y, p1.y, p2.y);

  // Clip agains screen boundries
  minX = std::max(minX, 0);
  minY = std::max(minY, 0);
  maxX = std::min(maxX, (int) screen_width - 1);
  maxY = std::min(maxY, (int) screen_height - 1);

  int bias0 = is_top_left(p1, p2) ? 0 : -1;
  int bias1 = is_top_left(p2, p0) ? 0 : -1;
  int bias2 = is_top_left(p0, p1) ? 0 : -1;

  // How to define point p in a way that makes sense?
  Point p;
  for (p.y = minY; p.y < maxY; p.y++)
  {
    for (p.x = minX; p.x < maxX; p.x++)
    {
      
      int w0 = orient2d(p1, p2, p) + bias0;
      int w1 = orient2d(p2, p0, p) + bias1;
      int w2 = orient2d(p0, p1, p) + bias2;
      if (w0 >= 0 && w1 >= 0 && w2 >= 0)
      {
        setImagePixel(screen_buffer, p, w0, w1, w2, p0.z, p1.z, p2.z, col);
      }
    }
  }
};

bool float_equal(float a, float b)
{
  return ((a - b) < 0.0001f && (a -b) > -0.0001f) ? true : false;
};

// Sanity check tests
void assert_matrix()
{
  Matrix4 m_identity;
  assert(float_equal(m_identity.m_data[0][0], 1));
  assert(float_equal(m_identity.m_data[1][1], 1));
  assert(float_equal(m_identity.m_data[2][2], 1));
  assert(float_equal(m_identity.m_data[3][3], 1));

  Matrix4 mul_a {
    5, 7, 9, 10,
    2, 3, 3, 8,
    8, 10, 2, 3,
    3, 3, 4, 8
  };
  Matrix4 mul_b {
    3, 10, 12, 18,
    12, 1, 4, 9,
    9, 10, 12, 2,
    3, 12, 4, 10
  };

  Matrix4 mul_res = GeometryCalc::matrix_matrix_mul(mul_a, mul_b);

  assert(float_equal(mul_res.m_data[0][0], 210.0));
  assert(float_equal(mul_res.m_data[0][1], 267.0));
  assert(float_equal(mul_res.m_data[0][2], 236.0));
  assert(float_equal(mul_res.m_data[0][3], 271.0));
  assert(float_equal(mul_res.m_data[1][1], 149.0));
  assert(float_equal(mul_res.m_data[2][2], 172.0));
  assert(float_equal(mul_res.m_data[3][0], 105.0));
  assert(float_equal(mul_res.m_data[3][3], 169.0));
  std::cout << "Asserting matrix: all good" << std::endl;
};

void assert_vectors()
{
    Vec4 v1{ 1.0, 1.0, 1.0, 1.0 };
    Vec4 v2{ 1.0, 1.0, 1.0, 1.0 };
    float res = GeometryCalc::vector_dot(v1, v2);
    assert(float_equal(res, 3.0));
};

