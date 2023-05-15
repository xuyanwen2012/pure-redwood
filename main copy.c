#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "float4.h"

float4 generate_random_float4(void) {
  float4 random_float4;
  random_float4.x = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.y = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.z = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.w = (float)rand() / RAND_MAX * 1000.0f;
  return random_float4;
}

// --------------- Utils ----------------------

void swap(float4* a, float4* b) {
  const float4 temp = *a;
  *a = *b;
  *b = temp;
}

int compare_dim(const float4 p1, const float4 p2, const int dim) {
  if (dim == 0) return p1.x - p2.x;
  if (dim == 1) return p1.y - p2.y;
  if (dim == 2) return p1.z - p2.z;
  return p1.w - p2.w;
}

int get_dim(const float4 p, const int dim) {
  if (dim == 0) return p.x;
  if (dim == 1) return p.y;
  if (dim == 2) return p.z;
  return p.w;
}

int partition(float4 arr[], const int low, const int high, const int dim) {
  int i = low - 1;

  for (int j = low; j <= high - 1; j++) {
    if (compare_dim(arr[j], arr[high], dim)) {
      i++;
      swap(&arr[i], &arr[j]);
    }
  }
  swap(&arr[i + 1], &arr[high]);
  return i + 1;
}

void nth_element(float4 arr[], const int low, const int high, const int n,
                 const int dim) {
  if (low < high) {
    const int pivot = partition(arr, low, high, dim);

    if (pivot == n) return;
    if (n < pivot)
      nth_element(arr, low, pivot - 1, n, dim);
    else
      nth_element(arr, pivot + 1, high, n, dim);
  }
}

void iota(int* array, int size, int start) {
  for (int i = 0; i < size; i++) {
    array[i] = start++;
  }
}

// -----------------  Tree --------------------------

enum {
  N = 1024,
  M = 8,
  MAX_NODES = 2048,
  MAX_STACK_SIZE = 128,
  LEAF_SIZE = 32,
};

typedef struct {
  int left;
  int right;
} Node;

float4 in_data[N];

int v_acc[MAX_NODES];
int next_node = 0;
Node children[MAX_NODES];

// int new_node(const float4 point) {
//   if (next_node < MAX_NODES) {
//     const int node_uid = next_node;
//     // node_data[node_uid] = point;
//     children[node_uid].left = -1;
//     children[node_uid].right = -1;
//     ++next_node;
//     return node_uid;
//   }

//   printf("node exceed");
//   return -1;
// }

int new_node(const int left, const int right) {
  if (next_node < MAX_NODES) {
    const int node_uid = next_node;
    children[node_uid].left = left;
    children[node_uid].right = right;
    ++next_node;
    return node_uid;
  }

  printf("node exceed");
  return -1;
}

int build_tree(const int left_idx, const int right_idx, const int depth) {
  int node_id = new_node(left_idx, right_idx);

  if (right_idx - left_idx <= LEAF_SIZE) {
  } else {
    
  }
}

void print_float4(const float4 p) {
  printf("(%f,\t%f,\t%f,\t%f)", p.x, p.y, p.w, p.z);
}

int main(void) {
  srand(114514);  // NOLINT(cert-msc51-cpp)

  for (int i = 0; i < N; i++) {
    in_data[i] = generate_random_float4();
  }

  int start = 0;
  int end = N;
  int mid = (start + end) / 2;

  int num_nodes_created = 1024;

  // // nth_element(dat1, start, end - 1, mid, 0);
  // for (int i = 0; i < N; ++i) {
  //   printf("%d:\t", i);
  //   print_float4(in_data[i]);
  //   printf("\n");
  // }

  for (int i = 0; i < num_nodes_created; ++i) {
    printf("%d:\t%d\n", i, v_acc[i]);
  }

  return 0;
}
