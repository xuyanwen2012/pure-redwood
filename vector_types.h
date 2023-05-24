#pragma once

#include <stdlib.h>

typedef struct float4 {
  float x, y, z, w;
} float4;

typedef struct float3 {
  float x, y, z;
} float3;

inline float3 make_float3(const float x, const float y, const float z) {
  float3 point;
  point.x = x;
  point.y = y;
  point.z = z;
  return point;
}

inline float4 make_float4(const float x, const float y, const float z,
                          const float w) {
  float4 point;
  point.x = x;
  point.y = y;
  point.z = z;
  point.w = w;
  return point;
}

inline float4 generate_random_float4(void) {
  float4 random_float4;
  random_float4.x = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.y = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.z = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.w = 1.0f;
  return random_float4;
}

inline void float4_add_assign(float4* a, const float4 b) {
  a->x += b.x;
  a->y += b.y;
  a->z += b.z;
  a->w += b.w;
}

inline float4 float4_divide(const float4* a, const float scalar) {
  return make_float4(a->x / scalar, a->y / scalar, a->z / scalar,
                     a->w / scalar);
}