#pragma once

#include <cstdint>
#include <cstring> // memcopy
//#include <cmath>

struct Vec4
{
  float x;
  float y;
  float z;
  float w;
};

struct Point
{
  //Vec4 in_space;
  float x;
  float y;
  float z;

  uint32_t r;
  uint32_t g;
  uint32_t b;
  uint32_t a;
};

struct Matrix4
{
  float m_data[4][4] {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  };
};

void to_identity_matrix(Matrix4 &m)
{
  m.m_data[0][0] = 1; m.m_data[0][1] = 0; m.m_data[0][2] = 0; m.m_data[0][3] = 0;
  m.m_data[1][0] = 0; m.m_data[1][1] = 1; m.m_data[1][2] = 0; m.m_data[1][3] = 0; 
  m.m_data[2][0] = 0; m.m_data[2][1] = 0; m.m_data[2][2] = 1; m.m_data[2][3] = 0; 
  m.m_data[3][0] = 0; m.m_data[3][1] = 0; m.m_data[3][2] = 0; m.m_data[3][3] = 1; 
};

struct Triangle
{
  Vec4 tri[3];
  float shadow_interpolated;
};

struct TriangleBuffer
{
  uint32_t max_num_tris;
  uint32_t size;
  Triangle *tris;
};

struct ScreenBuffer
{
  unsigned int r, g, b;
  float z;
};

void triangle_buffer_init(TriangleBuffer &t_buffer, uint32_t size, uint32_t max_num_tris)
{
  t_buffer.size = size;
  t_buffer.max_num_tris = max_num_tris;
  t_buffer.tris = (Triangle *) malloc(t_buffer.max_num_tris * sizeof(Triangle));
  if (t_buffer.tris == NULL)
    exit(-1);
  // t_buffer.tris = (Triangle *) malloc(t_buffer.max_num_tris * sizeof(Triangle));
};

void triangle_buffer_copy(TriangleBuffer *from_buffer, TriangleBuffer *to_buffer)
{
  to_buffer->max_num_tris = from_buffer->max_num_tris;
  to_buffer->size = from_buffer->size;
  memcpy (to_buffer->tris, from_buffer->tris, from_buffer->max_num_tris * sizeof(Triangle));
}

void triangle_buffer_delete(TriangleBuffer &t_buffer)
{
  free(t_buffer.tris);
  std::cout << "Free world_buffer" << std::endl;  
};

void triangle_buffer_add_triangle(TriangleBuffer &t_buffer, Triangle t)
{
  t_buffer.size++;
  if (t_buffer.size <= t_buffer.max_num_tris)
  {
    t_buffer.tris[t_buffer.size - 1] = t;
  }
}

struct GeometryCalc
{
  static Matrix4 matrix_matrix_mul(const Matrix4 &, const Matrix4 &);
  static Matrix4 matrix_matrix_add(const Matrix4 &, const Matrix4 &);
  static void matrix_row_to_column(Matrix4 &);

  static Vec4 vec_matrix_mul(const Vec4 &, const Matrix4 &);

  static Vec4 triangle_normal(const Triangle &a);

  // Vector math
  static Vec4 vector_add(const Vec4&, const Vec4&);
  static Vec4 vector_sub(const Vec4&, const Vec4&);
  static float vector_dot(const Vec4&, const Vec4&);
  static Vec4 vector_cross(const Vec4 &, const Vec4 &);
  static void normalize(Vec4 &);
};

// a is row/col matrix, b is col/row
Matrix4 GeometryCalc::matrix_matrix_mul(const Matrix4 &a, const Matrix4 &b)
{
  Matrix4 ret;
  
  for (int i = 0; i < 4; i++)
  {
      for (int j = 0; j < 4; j++)
      {
          ret.m_data[i*3][j] = 0;
          for (int k = 0; k < 4; k++)
          {
              ret.m_data[i * 3][j] += a.m_data[i * 3][k] * b.m_data[i * 3][k];
          }
      }
  }
  return ret;
};

Vec4 GeometryCalc::vec_matrix_mul(const Vec4 &v, const Matrix4 &m)
{
  Vec4 new_vec;
  new_vec.x = v.x * m.m_data[0][0] + v.y * m.m_data[1][0] + v.z * m.m_data[2][0] + v.w * m.m_data[3][0];
  new_vec.y = v.x * m.m_data[0][1] + v.y * m.m_data[1][1] + v.z * m.m_data[2][1] + v.w * m.m_data[3][1];
  new_vec.z = v.x * m.m_data[0][2] + v.y * m.m_data[1][2] + v.z * m.m_data[2][2] + v.w * m.m_data[3][2];
  new_vec.w = v.x * m.m_data[0][3] + v.y * m.m_data[1][3] + v.z * m.m_data[2][3] + v.w * m.m_data[3][3];
  return new_vec;
};

Matrix4 GeometryCalc::matrix_matrix_add(const Matrix4 &a, const Matrix4 &b)
{
  Matrix4 ret;
  for (uint32_t row = 0; row < 4; row++)
    for (uint32_t col = 0; col < 4; col++)
      ret.m_data[row][col] = a.m_data[row][col] + b.m_data[row][col];
  return ret;
}

void GeometryCalc::matrix_row_to_column(Matrix4 &matrix)
{
  float temp;
  for(size_t row = 0; row < 4; row++) {
    for(size_t col = 0; col < 4; col++) {
      temp = matrix.m_data[col][row];
      matrix.m_data[row][col] = temp;
    }
  }
}

Vec4 GeometryCalc::triangle_normal(const Triangle &t)
{
    Vec4 a = GeometryCalc::vector_sub(t.tri[1], t.tri[0]);
    Vec4 b = GeometryCalc::vector_sub(t.tri[2], t.tri[0]);
  
    GeometryCalc::normalize(a);
    GeometryCalc::normalize(b);

    return GeometryCalc::vector_cross(a, b);
}

// VECTOR 

Vec4 GeometryCalc::vector_add(const Vec4 &a, const Vec4 &b)
{
  Vec4 res;
  res.x = a.x + b.x;
  res.y = a.y + b.y;
  res.z = a.z + b.z;
  res.w = a.w; // wtf to do with w?
  return res;
}

Vec4 GeometryCalc::vector_sub(const Vec4 &a, const Vec4 &b)
{
  Vec4 res;
  res.x = a.x - b.x;
  res.y = a.y - b.y;
  res.z = a.z - b.z;
  res.w = a.w; // wtf to do with w?
  return res;
}

float GeometryCalc::vector_dot(const Vec4 &a, const Vec4 &b)
{
  return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

Vec4 GeometryCalc::vector_cross(const Vec4 &a, const Vec4 &b)
{
  // a2b3−a3b2,a3b1−a1b3,a1b2−a2b1
  Vec4 res;
  res.x = (a.y * b.z) - (a.z * b.y);
  res.y = (a.z * b.x) - (a.x * b.z);
  res.z = (a.x * b.y) - (a.y * b.x);
  res.w = a.w; // wtf to do with w?
  return res;
}

void GeometryCalc::normalize(Vec4 &vec)
{
  float magnitude = sqrtf((vec.x * vec.x) + (vec.y * vec.y) + (vec.z * vec.z));
  vec.x = vec.x / magnitude;
  vec.y = vec.y / magnitude;
  vec.z = vec.z / magnitude;
}