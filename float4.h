#ifndef PURE_REDWOOD_FLOAT4_H
#define PURE_REDWOOD_FLOAT4_H

struct float4 {
  float x, y, z, w;
};

typedef struct float4 float4;

inline float4 make_float4(const float x, const float y, const float z,
                          const float w) {
  float4 result;
  result.x = x;
  result.y = y;
  result.z = z;
  result.w = w;
  return result;
}

inline float4 add_f4(const float4 a, const float4 b) {
  float4 result;
  result.x = a.x + b.x;
  result.y = a.y + b.y;
  result.z = a.z + b.z;
  result.w = a.w + b.w;
  return result;
}

inline void add_assign_f4(float4* const a, const float4 b) {
  a->x += b.x;
  a->y += b.y;
  a->z += b.z;
  a->w += b.w;
}

inline float4 subtract_f4(const float4 a, const float4 b) {
  float4 result;
  result.x = a.x - b.x;
  result.y = a.y - b.y;
  result.z = a.z - b.z;
  result.w = a.w - b.w;
  return result;
}

inline void subtract_assign_f4(float4* const a, const float4 b) {
  a->x -= b.x;
  a->y -= b.y;
  a->z -= b.z;
  a->w -= b.w;
}

inline float4 multiply_f4(const float4 a, const float4 b) {
  float4 result;
  result.x = a.x * b.x;
  result.y = a.y * b.y;
  result.z = a.z * b.z;
  result.w = a.w * b.w;
  return result;
}

inline void multiply_assign_f4(float4* const a, const float4 b) {
  a->x *= b.x;
  a->y *= b.y;
  a->z *= b.z;
  a->w *= b.w;
}

inline float4 divide_f4(const float4 a, const float4 b) {
  float4 result;
  result.x = a.x / b.x;
  result.y = a.y / b.y;
  result.z = a.z / b.z;
  result.w = a.w / b.w;
  return result;
}

inline void divide_assign_f4(float4* const a, const float4 b) {
  a->x /= b.x;
  a->y /= b.y;
  a->z /= b.z;
  a->w /= b.w;
}

#endif  // PURE__REDWOOD_FLOAT4_H
