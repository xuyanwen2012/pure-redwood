#include <float.h>
#include <math.h>
#include <stddef.h>
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

// -----------------  Tree --------------------------

enum {
  N = 1024,
  M = 8,
  MAX_NODES = 2048,
  MAX_STACK_SIZE = 128,
};

typedef struct BranchNode {
  int left;
  int right;
  float4 data;
} BranchNode;

typedef struct LeafNode {
  int head_addr;
  int count;
} LeafNode;

typedef enum { BRANCH, LEAF } NodeType;

typedef struct TreeNode {
  NodeType type;

  union {
    BranchNode branch;
    LeafNode leaf;
  };
} TreeNode;

// ----------------------
float4 in_data[N];

int next_node = 0;
TreeNode nodes[MAX_NODES];

// ----------------------

int new_branch_node(void) {
  const int cur = next_node;
  nodes[cur].type = BRANCH;
  nodes[cur].branch.left = -1;
  nodes[cur].branch.right = -1;
  ++next_node;
  return cur;
}

int new_leaf_node(void) {
  const int cur = next_node;
  nodes[cur].type = LEAF;
  nodes[cur].leaf.head_addr = -1;
  nodes[cur].leaf.count = 0;
  ++next_node;
  return cur;
}

void print_float4(const float4 p) {
  printf("(%f,\t%f,\t%f,\t%f)", p.x, p.y, p.w, p.z);
}

int build_tree(float4* data, const int low, const int high, const int depth,
               const int max_leaf_size) {
  int node_id;

  const int count = high - low;
  if (count <= max_leaf_size) {
    node_id = new_leaf_node();

    nodes[node_id].leaf.count = count;
    nodes[node_id].leaf.head_addr = low;

  } else {
    node_id = new_branch_node();

    const int axis = depth % 4;
    const int mid = (low + high) / 2;

    nth_element(data, low, high, mid, axis);

    nodes[node_id].branch.data = data[mid];
    nodes[node_id].branch.left =
        build_tree(data, low, mid, depth + 1, max_leaf_size);
    nodes[node_id].branch.right =
        build_tree(data, mid + 1, high, depth + 1, max_leaf_size);
  }

  return node_id;
}

void dfs(int cur, const int depth) {
  if (cur == -1) {
    return;
  } else if (nodes[cur].type == LEAF) {
    printf("(l) head:%d\tcount:%d\n", nodes[cur].leaf.head_addr,
           nodes[cur].leaf.count);
  } else if (nodes[cur].type == BRANCH) {
    dfs(nodes[cur].branch.left, depth + 1);
    dfs(nodes[cur].branch.right, depth + 1);

    const int axis = depth % 4;
    printf("(b) axis:%d\tpoint:", axis);
    print_float4(nodes[cur].branch.data);
    printf("\n");
  }
}

int main(void) {
  srand(114514);  // NOLINT(cert-msc51-cpp)

  for (int i = 0; i < N; i++) {
    in_data[i] = generate_random_float4();
  }

  // we want to make this a dynamic parameter
  const int max_leaf_size = 32;
  const int root = build_tree(in_data, 0, N, 0, max_leaf_size);

  dfs(root, 0);

  printf("Nodes used: %d\n", next_node);

  return 0;
}
