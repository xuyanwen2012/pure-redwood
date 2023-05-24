#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "vector_types.h"

// ---------------------------------------------------------------------------
//  App Constants
// ---------------------------------------------------------------------------

enum {
  N = 10240,
  M = 8,
  MAX_NODES = 1024 * 5,
  MAX_STACK_SIZE = 64,

  TREE_BUILD_LEAF_SIZE = 16,

  // how many elements can be processed by the FPGA
  DUET_LEAF_SIZE = 32,
  NUM_EXECUTORS = 2,
};

#define SOFTENING 1e-9f
#define THETA 0.2f

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

    int num_valid = 0;
    for (int i = 0; i < 8; ++i) {
      const int child = nodes[cur].children[i];

      if (!isnan(nodes[child].center_of_mass.x)) {
        float4_add_assign(&com, nodes[child].center_of_mass);
        ++num_valid;
      }
    }

    nodes[cur].center_of_mass = float4_divide(&com, num_valid);
  }
}

float compute_theta_value(const int node_id, const float3 q) {
  const float4 com = nodes[node_id].center_of_mass;
  float norm_sqr = 1e-9f;
  const float dx = com.x - q.x;
  const float dy = com.y - q.y;
  const float dz = com.z - q.z;
  norm_sqr += dx * dx + dy * dy + dz * dz;
  const float norm = sqrtf(norm_sqr);
  const float geo_size = boxes[node_id].x_max - boxes[node_id].x_min;
  return geo_size / norm;
}

float gravity_func_4f(const float4 p, const float3 q) {
  const float dx = p.x - q.x;
  const float dy = p.y - q.y;
  const float dz = p.z - q.z;
  const float dist_sqr = dx * dx + dy * dy + dz * dz + SOFTENING;
  const float inv_dist = 1.0f / sqrtf(dist_sqr);
  const float inv_dist3 = inv_dist * inv_dist * inv_dist;
  const float with_mass = inv_dist3 * p.w;
  return dx * with_mass + dy * with_mass + dz * with_mass;
}

int stats_num_branch_visited = 0;
int stats_num_leaf_visited = 0;
int stats_num_elements = 0;

void compute_gravity_at(const float4* in, const int cur, const float3 q,
                        float* sum) {
  if (compute_theta_value(cur, q) < THETA) {
    *sum += gravity_func_4f(nodes[cur].center_of_mass, q);
    ++stats_num_elements;
  } else if (nodes[cur].is_leaf) {
    ++stats_num_leaf_visited;
    for (int i = ranges[cur].low; i < ranges[cur].high; ++i) {
      *sum += gravity_func_4f(in[i], q);
      ++stats_num_elements;
    }
  } else {
    ++stats_num_branch_visited;
    for (int i = 0; i < 8; ++i) {
      compute_gravity_at(in, nodes[cur].children[i], q, sum);
    }
  }
}

int main(void) {
  for (int i = 0; i < N; ++i) in_data[i] = generate_random_float4();

  const int root =
      build_tree(in_data, 0, N, make_box_3d(0, 0, 0, 1000, 1000, 1000), 0);

  compute_center_of_masses(in_data, root, 0);

  printf("%f\n", nodes[root].center_of_mass.x);
  printf("%f\n", nodes[root].center_of_mass.y);
  printf("%f\n", nodes[root].center_of_mass.z);
  printf("%f\n", nodes[root].center_of_mass.w);

  float sum = 0.0f;
  const float3 q = make_float3(500, 500, 500);

  compute_gravity_at(in_data, root, q, &sum);

  printf("sum: %f\n", sum);

  float gt = 0.0f;
  for (int i = 0; i < N; ++i) {
    gt += gravity_func_4f(in_data[i], q);
  }

  printf("gt: %f\n", gt);

  printf("stats_num_branch_visited: %d\n", stats_num_branch_visited);
  printf("stats_num_leaf_visited: %d\n", stats_num_leaf_visited);
  printf("stats_num_elements: %d\n", stats_num_elements);

  return EXIT_SUCCESS;
}
