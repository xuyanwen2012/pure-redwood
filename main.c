#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "float4.h"

#define REDWOOD_DEBUG

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

// // We want user write this
// typedef struct {
//   int left;
//   int right;
//   float4 data;
// } Node;

// We convert to this

// typedef struct {
//   int left;
//   int right;
// } BranchNode;

typedef struct BranchNode {
  struct TreeNode* left;
  struct TreeNode* right;
  int axis;
  float4 data;
} BranchNode;

typedef struct LeafNode {
  float4* head_addr;
  int count;
} LeafNode;

typedef enum { BRANCH, LEAF } NodeType;

typedef struct TreeNode {
  NodeType type;

  // DEBUG ONLY
#ifdef REDWOOD_DEBUG
  int uid;
#endif

  union {
    BranchNode branch;
    LeafNode leaf;
  };
} TreeNode;

float4 in_data[N];

int v_acc[MAX_NODES];
int next_node = 0;

TreeNode* new_branch_node(void) {
  TreeNode* node = malloc(sizeof(TreeNode));
  node->type = BRANCH;
#ifdef REDWOOD_DEBUG
  node->uid = next_node++;
#endif
  node->branch.left = NULL;
  node->branch.right = NULL;
  node->branch.axis = -1;
  node->branch.data.x = -1.0f;
  node->branch.data.y = -1.0f;
  node->branch.data.z = -1.0f;
  node->branch.data.w = -1.0f;
  return node;
}

TreeNode* new_leaf_node(void) {
  TreeNode* node = malloc(sizeof(TreeNode));
  node->type = LEAF;
#ifdef REDWOOD_DEBUG
  node->uid = next_node++;
#endif
  node->leaf.head_addr = NULL;
  node->leaf.count = 0;
  return node;
}

int insert_to_leaf(LeafNode* node, const float4 p) {
  const int cur = node->count;

  if (cur == 0) {
    node->head_addr = malloc(sizeof(float4) * LEAF_SIZE);
  }

  if (cur < LEAF_SIZE) {
    node->head_addr[cur] = p;
    ++node->count;
    return 0;
  }
  // Need to split
  return -1;
}

void print_float4(const float4 p) {
  printf("(%f,\t%f,\t%f,\t%f)", p.x, p.y, p.w, p.z);
}

void print_leaf(const LeafNode* node) {
  const int num = node->count;
  printf("num: %d\n", num);
  for (int i = 0; i < num; ++i) {
    printf("\t%d:\t", i);
    print_float4(node->head_addr[i]);
    printf("\n");
  }
  printf("\n");
}

TreeNode* build_tree(float4* data, const int low, const int high,
                     const int depth) {
  TreeNode* node = NULL;

  const int count = high - low;
  if (count <= LEAF_SIZE) {
    node = new_leaf_node();

    node->leaf.count = count;
    node->leaf.head_addr = &data[low];

  } else {
    node = new_branch_node();

    const int axis = depth % 4;
    const int mid = (low + high) / 2;

    nth_element(data, low, high, mid, axis);

    node->branch.data = data[mid];
    node->branch.axis = axis;
    node->branch.left = build_tree(data, low, mid, depth + 1);
    node->branch.right = build_tree(data, mid + 1, high, depth + 1);
  }

  return node;
}

void dfs(TreeNode* node) {
  if (node == NULL) {
    return;
  } else if (node->type == LEAF) {
    float4* my = node->leaf.head_addr;
    ptrdiff_t diff = my - &in_data[0];

    printf("(l) head:%p (%td)\tcount:%d\n", (void*)node->leaf.head_addr, diff,
           node->leaf.count);
  } else if (node->type == BRANCH) {
    printf("(b) axis:%d\tpoint:", node->branch.axis);
    print_float4(node->branch.data);
    printf("\n");

    dfs(node->branch.left);
    dfs(node->branch.right);
  }
}

int main(void) {
  srand(114514);  // NOLINT(cert-msc51-cpp)

  for (int i = 0; i < N; i++) {
    in_data[i] = generate_random_float4();
  }

  TreeNode* root = build_tree(in_data, 0, N, 0);

  for (int i = 0; i < N; ++i) {
    printf("%d:\t", i);
    print_float4(in_data[i]);
    printf("\n");
  }

  dfs(root);

  return 0;
}
