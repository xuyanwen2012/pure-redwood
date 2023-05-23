#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// ---------------------------------------------------------------------------
//  App Constants
// ---------------------------------------------------------------------------

#ifdef REDWOOD_DEBUG
#define DEBUG_PRINT_DASH(depth) \
  for (int j = 0; j < depth; ++j) putchar('-');
#define DEBUG_PRINT(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT_DASH(depth) \
  do {                          \
  } while (0)
#define DEBUG_PRINT(fmt, ...) \
  do {                        \
  } while (0)
#endif

enum {
  N = 10240,
  M = 8,
  MAX_NODES = 2048,
  MAX_STACK_SIZE = 64,

  TREE_BUILD_LEAF_SIZE = 32,

  // how many elements can be processed by the FPGA
  DUET_LEAF_SIZE = 32,
  NUM_EXECUTORS = 2,
};

// ---------------------------------------------------------------------------
//  Utils
// ---------------------------------------------------------------------------

typedef struct float4 {
  float x, y, z, w;
} float4;

typedef struct float3 {
  float x, y, z;
} float3;

float3 make_float3(const float x, const float y, const float z) {
  float3 point;
  point.x = x;
  point.y = y;
  point.z = z;
  return point;
}

float4 make_float4(const float x, const float y, const float z, const float w) {
  float4 point;
  point.x = x;
  point.y = y;
  point.z = z;
  point.w = w;
  return point;
}

float4 generate_random_float4(void) {
  float4 random_float4;
  random_float4.x = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.y = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.z = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.w = 1.0f;
  return random_float4;
}

void float4_add_assign(float4* a, const float4 b) {
  a->x += b.x;
  a->y += b.y;
  a->z += b.z;
  a->w += b.w;
}

float4 float4_divide(const float4* a, const float scalar) {
  return make_float4(a->x / scalar, a->y / scalar, a->z / scalar,
                     a->w / scalar);
}

// ---------------------------------------------------------------------------
//  Bounding Box 3D
// ---------------------------------------------------------------------------

typedef struct BoundingBox3D {
  float x_min, y_min, z_min, x_max, y_max, z_max;
} BoundingBox3D;

BoundingBox3D make_box_3d(const float xmin, const float ymin, const float zmin,
                          const float xmax, const float ymax,
                          const float zmax) {
  BoundingBox3D box;
  box.x_min = xmin;
  box.y_min = ymin;
  box.z_min = zmin;
  box.x_max = xmax;
  box.y_max = ymax;
  box.z_max = zmax;
  return box;
}

// (0,0) - (1, 0) - ... (row, 0)
// (0,1) - (1, 1) - ... (row, 1)
// ...
// (0,c) - (1, c) - ... (row, col)
//

// (x_min, y_min) - (x_mid, y_min) - (x_max, y_min)
// (x_min, y_mid) - (x_mid, y_mid) - (x_max, y_mid)
// (x_min, y_max) - (x_mid, y_max) - (x_max, y_max)
void split_box(const BoundingBox3D* box, BoundingBox3D* sub_boxes) {
  const float x_min = box->x_min;
  const float y_min = box->y_min;
  const float z_min = box->z_min;
  const float x_max = box->x_max;
  const float y_max = box->y_max;
  const float z_max = box->z_max;
  const float x_mid = (x_min + x_max) / 2.0f;
  const float y_mid = (y_min + y_max) / 2.0f;
  const float z_mid = (z_min + z_max) / 2.0f;
  sub_boxes[0] = make_box_3d(x_min, y_min, z_min, x_mid, y_mid, z_mid);
  sub_boxes[1] = make_box_3d(x_mid, y_min, z_min, x_max, y_mid, z_mid);
  sub_boxes[2] = make_box_3d(x_min, y_mid, z_min, x_mid, y_max, z_mid);
  sub_boxes[3] = make_box_3d(x_mid, y_mid, z_min, x_max, y_max, z_mid);
  sub_boxes[4] = make_box_3d(x_min, y_min, z_mid, x_mid, y_mid, z_max);
  sub_boxes[5] = make_box_3d(x_mid, y_min, z_mid, x_max, y_mid, z_max);
  sub_boxes[6] = make_box_3d(x_min, y_mid, z_mid, x_mid, y_max, z_max);
  sub_boxes[7] = make_box_3d(x_mid, y_mid, z_mid, x_max, y_max, z_max);
}

// int determine_quadrant(const BoundingBox3D* box, const float x, const float
// y,
//                        const float z) {
//   const float x_mid = (box->x_min + box->x_max) / 2.0f;
//   const float y_mid = (box->y_min + box->y_max) / 2.0f;
//   const float z_mid = (box->z_min + box->z_max) / 2.0f;

//   int quadrant;
//   if (x < x_mid) {
//     if (y < y_mid) {
//       if (z < z_mid) {
//         quadrant = 0;
//       } else {
//         quadrant = 1;
//       }
//     } else if (z < z_mid) {
//       quadrant = 2;
//     } else {
//       quadrant = 3;
//     }
//   } else if (y < y_mid) {
//     if (z < z_mid) {
//       quadrant = 4;
//     } else {
//       quadrant = 5;
//     }
//   } else if (z < z_mid) {
//     quadrant = 6;
//   } else {
//     quadrant = 7;
//   }
//   return quadrant;
// }

int determine_quadrant(const BoundingBox3D* box, const float x, const float y,
                       const float z) {
  const float x_mid = (box->x_min + box->x_max) / 2.0f;
  const float y_mid = (box->y_min + box->y_max) / 2.0f;
  const float z_mid = (box->z_min + box->z_max) / 2.0f;

  int quadrant;
  // Determine the quadrant based on the point's position relative
  // to the box's midpoint
  if (x < x_mid) {
    if (y < y_mid) {
      if (z < z_mid) {
        quadrant = 0;  // Bottom-back-left quadrant
      } else {
        quadrant = 4;  // Top-back-left quadrant
      }
    } else {
      if (z < z_mid) {
        quadrant = 2;  // Bottom-front-left quadrant
      } else {
        quadrant = 6;  // Top-front-left quadrant
      }
    }
  } else {
    if (y < y_mid) {
      if (z < z_mid) {
        quadrant = 1;  // Bottom-back-right quadrant
      } else {
        quadrant = 5;  // Top-back-right quadrant
      }
    } else {
      if (z < z_mid) {
        quadrant = 3;  // Bottom-front-right quadrant
      } else {
        quadrant = 7;  // Top-front-right quadrant
      }
    }
  }

  return quadrant;
}

// ---------------------------------------------------------------------------
//  Tree Node
// ---------------------------------------------------------------------------

typedef struct Node {
  bool is_leaf;
  float4 center_of_mass;
  int children[8];
} Node;

typedef struct Range {
  int low;
  int high;
} Range;

// All the data are here
BoundingBox3D boxes[MAX_NODES];
Node nodes[MAX_NODES];
Range ranges[MAX_NODES];
float4 in_data[N];

int next_node = 0;
int get_next_node_id(void) { return next_node++; }

typedef struct {
  float x, y, z;
  BoundingBox3D* box_ptr;
} Args;

const BoundingBox3D* tmp_box_ptr;

int compare(const void* a, const void* b) {
  const float4 arg1 = *(const float4*)a;
  const float4 arg2 = *(const float4*)b;
  const int group_id_a =
      determine_quadrant(tmp_box_ptr, arg1.x, arg1.y, arg1.z);
  const int group_id_b =
      determine_quadrant(tmp_box_ptr, arg2.x, arg2.y, arg2.z);
  if (group_id_a < group_id_b) return -1;
  if (group_id_a > group_id_b) return 1;
  return 0;
}

int build_tree(float4* in, const int low, const int high,
               const BoundingBox3D box, const int depth) {
  const int cur = get_next_node_id();
  boxes[cur] = box;

  const int n = high - low;
  if (n <= TREE_BUILD_LEAF_SIZE) {
    nodes[cur].is_leaf = true;
    ranges[cur].low = low;
    ranges[cur].high = high;
  } else {
    nodes[cur].is_leaf = false;
    BoundingBox3D sub_boxes[8];
    split_box(&box, &sub_boxes[0]);

    tmp_box_ptr = &box;
    qsort(in + low, n, sizeof(float4), compare);

    volatile int count[8];
    for (int i = 0; i < 8; ++i) count[i] = 0;

    for (int i = low; i < high; ++i) {
      const float x = in[i].x;
      const float y = in[i].y;
      const float z = in[i].z;
      const int quadrant = determine_quadrant(&box, x, y, z);
      ++count[quadrant];
    }

    int next_low = low;
    for (int i = 0; i < 8; ++i) {
      const int next_high = next_low + count[i];

      DEBUG_PRINT_DASH(depth);
      DEBUG_PRINT("[%d, %d) - %d\n", next_low, next_high, count[i]);

      nodes[cur].children[i] =
          build_tree(in, next_low, next_high, sub_boxes[i], depth + 1);

      next_low = next_high;
    }

    ranges[cur].low = ranges[nodes[cur].children[0]].low;
    ranges[cur].high = ranges[nodes[cur].children[7]].high;
  }

  return cur;
}

void compute_center_of_masses(float4* in, const int cur, const int depth) {
  float4 com = make_float4(0.0f, 0.0f, 0.0f, 0.0f);

  if (nodes[cur].is_leaf) {
    for (int i = ranges[cur].low; i < ranges[cur].high; ++i) {
      float4_add_assign(&com, in[i]);
    }
    const int count = ranges[cur].high - ranges[cur].low;
    nodes[cur].center_of_mass = float4_divide(&com, count);

  } else {
    for (int i = 0; i < 8; ++i) {
      const int child = nodes[cur].children[i];
      compute_center_of_masses(in, child, depth + 1);
    }

    for (int i = 0; i < 8; ++i) {
      const int child = nodes[cur].children[i];
      float4_add_assign(&com, nodes[child].center_of_mass);
    }

    nodes[cur].center_of_mass = float4_divide(&com, 8);
  }
}

int main(void) {
  for (int i = 0; i < N; ++i) in_data[i] = generate_random_float4();

  const int root =
      build_tree(in_data, 0, N, make_box_3d(0, 0, 0, 1000, 1000, 1000), 0);

  compute_center_of_masses(in_data, root, 0);

  // printf("%f\n", nodes[0].center_of_mass.w);

  return EXIT_SUCCESS;
}
