#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "float4.h"

// ---------------------------------------------------------------------------
//  Utils
// ---------------------------------------------------------------------------

float4 generate_random_float4(void) {
  float4 random_float4;
  random_float4.x = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.y = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.z = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.w = (float)rand() / RAND_MAX * 1000.0f;
  return random_float4;
}

enum {
  N = 1024 * 10,
  M = 8,
  MAX_NODES = 2048,
  MAX_STACK_SIZE = 128,
};

typedef struct BranchNode {
  struct TreeNode* left_child;
  struct TreeNode* right_child;
  float4 point;
} BranchNode;

typedef struct LeafNode {
  float4* head_addr;
  int count;
} LeafNode;

typedef enum { BRANCH, LEAF } NodeType;

typedef enum { LEFT = 0, RIGHT = 1 } Direction;

typedef struct TreeNode {
  NodeType type;
  int uid;

  union {
    BranchNode branch;
    LeafNode leaf;
  };
} TreeNode;

int compare_dim(const float4 p1, const float4 p2, const int dim) {
  // Should be optimized away, double check
  if (dim == 0) return p1.x < p2.x;
  if (dim == 1) return p1.y < p2.y;
  if (dim == 2) return p1.z < p2.z;
  return p1.w < p2.w;
}

float get_dim(const int dim, const float4 p) {
  // Should be optimized away, double check
  if (dim == 0) return p.x;
  if (dim == 1) return p.y;
  if (dim == 2) return p.z;
  return p.w;
}

// ---------------------------------------------------------------------------
//  Kernel Functions
// ---------------------------------------------------------------------------

inline float kernel_func_4f(const float4 p, const float4 q) {
  const float dx = p.x - q.x;
  const float dy = p.y - q.y;
  const float dz = p.z - q.z;
  const float dw = p.w - q.w;
  return sqrtf(dx * dx + dy * dy + dz * dz + dw * dw);
}

inline float kernel_func_1f(const float p, const float q) {
  const float dx = p - q;
  return sqrtf(dx * dx);
}

inline void reduce_element(const float4 p, const float4 q, float* my_min) {
  const float dist = kernel_func_4f(p, q);
  *my_min = fminf(*my_min, dist);
}

float ground_truth(const float4* data, const float4 q) {
  float gt = FLT_MAX;
  for (int i = 0; i < N; ++i) {
    const float dist = kernel_func_4f(data[i], q);
    gt = fminf(gt, dist);
  }
  return gt;
}

// ---------------------------------------------------------------------------
//  C++ STD Algorithm
// ---------------------------------------------------------------------------

typedef float4* RanIt;
typedef float4* BidIt;

float4* prev_iter(float4* first) { return --first; }

float4* next_iter(float4* first) { return ++first; }

void iter_swap(float4* a, float4* b) {
  const float4 temp = *a;
  *a = *b;
  *b = temp;
}

void med_3_unchecked(const RanIt first, const RanIt mid, const RanIt last,
                     const int axis) {
  // sort median of three elements to middle
  if (compare_dim(*mid, *first, axis)) {
    iter_swap(mid, first);
  }

  if (compare_dim(*last, *mid, axis)) {
    iter_swap(last, mid);

    if (compare_dim(*mid, *first, axis)) {
      iter_swap(mid, first);
    }
  }
}

void guess_median_unchecked(const RanIt first, const RanIt mid,
                            const RanIt last, const int axis) {
  // sort median element to middle
  const ptrdiff_t count = last - first;

  if (40 < count) {
    // Tukey's ninther
    const ptrdiff_t step = (count + 1) >> 3;
    // made inclusive in caller
    const ptrdiff_t two_step = step << 1;
    med_3_unchecked(first, first + step, first + two_step, axis);
    med_3_unchecked(mid - step, mid, mid + step, axis);
    med_3_unchecked(last - two_step, last - step, last, axis);
    med_3_unchecked(first + step, mid, last - step, axis);
  } else {
    med_3_unchecked(first, mid, last, axis);
  }
}

typedef struct Pair {
  float4* first;
  float4* second;
} Pair;

Pair partition_by_median_guess_unchecked(const RanIt first, const RanIt last,
                                         const int axis) {
  // partition [_First, _Last)
  const RanIt mid = first + ((last - first) >> 1);  // shift for codegen
  guess_median_unchecked(first, mid, prev_iter(last), axis);
  RanIt pfirst = mid;
  RanIt plast = next_iter(pfirst);

  while (first < pfirst && !compare_dim(*prev_iter(pfirst), *pfirst, axis) &&
         !compare_dim(*pfirst, *prev_iter(pfirst), axis)) {
    --pfirst;
  }

  while (plast < last && !compare_dim(*plast, *pfirst, axis) &&
         !compare_dim(*pfirst, *plast, axis)) {
    ++plast;
  }

  RanIt gfirst = plast;
  RanIt glast = pfirst;

  for (;;) {
    // partition
    for (; gfirst < last; ++gfirst) {
      if (compare_dim(*pfirst, *gfirst, axis)) {
        continue;
      }
      if (compare_dim(*gfirst, *pfirst, axis)) {
        break;
      }
      if (plast != gfirst) {
        iter_swap(plast, gfirst);
        ++plast;
      } else {
        ++plast;
      }
    }

    for (; first < glast; --glast) {
      if (compare_dim(*prev_iter(glast), *pfirst, axis)) {
        continue;
      }
      if (compare_dim(*pfirst, *prev_iter(glast), axis)) {
        break;
      }
      if (--pfirst != prev_iter(glast)) {
        iter_swap(pfirst, prev_iter(glast));
      }
    }

    if (glast == first && gfirst == last) {
      Pair p;
      p.first = pfirst;
      p.second = plast;
      return p;
    }

    if (glast == first) {
      // no room at bottom, rotate pivot upward
      if (plast != gfirst) {
        iter_swap(pfirst, plast);
      }

      ++plast;
      iter_swap(pfirst, gfirst);
      ++pfirst;
      ++gfirst;
    } else if (gfirst == last) {
      // no room at top, rotate pivot downward
      if (--glast != --pfirst) {
        iter_swap(glast, pfirst);
      }

      iter_swap(pfirst, --plast);
    } else {
      iter_swap(gfirst, --glast);
      ++gfirst;
    }
  }
}

BidIt move_backward_unchecked(const BidIt first, BidIt last, BidIt dest) {
  while (first != last) {
    *--dest = *--last;
  }
  return dest;
}

BidIt insertion_sort_unchecked(const BidIt first, const BidIt last,
                               const int axis) {
  // insertion sort [_First, _Last)
  if (first != last) {
    for (BidIt mid = first; ++mid != last;) {
      // order next element
      BidIt hole = mid;
      const float4 val = *mid;

      if (compare_dim(val, *first, axis)) {
        // found new earliest element, move to front
        move_backward_unchecked(first, mid, ++hole);
        *first = val;
      } else {
        // look for insertion point after first
        for (BidIt prev = hole; compare_dim(val, *--prev, axis); hole = prev) {
          *hole = *prev;  // move hole down
        }

        *hole = val;  // insert element in hole
      }
    }
  }

  return last;
}

void nth_element(const RanIt first, const RanIt nth, const RanIt last,
                 const int axis) {
  // order Nth element
  RanIt u_first = first;
  const RanIt u_nth = nth;
  RanIt u_last = last;
  if (u_nth == u_last) {
    return;  // nothing to do
  }

  // ISORT_MAX = 32
  while (32 < u_last - u_first) {
    // divide and conquer, ordering
    // partition containing Nth
    const Pair u_mid =
        partition_by_median_guess_unchecked(u_first, u_last, axis);

    if (u_mid.second <= u_nth) {
      u_first = u_mid.second;
    } else if (u_mid.first <= u_nth) {
      return;  // _Nth is in the subrange of elements equal to the pivot; done
    } else {
      u_last = u_mid.first;
    }
  }

  insertion_sort_unchecked(u_first, u_last, axis);  // sort any remainder
}

// ---------------------------------------------------------------------------
//  Tree data structure
// ---------------------------------------------------------------------------

int next_node = 0;

TreeNode* BuildTree(float4* data, const int low, const int high,
                    const int leaf_size, const int depth) {
  TreeNode* node = malloc(sizeof(TreeNode));
  node->uid = next_node++;

  const int n = high - low;
  if (n <= leaf_size) {
    node->type = LEAF;

    node->leaf.head_addr = (float4*)malloc(sizeof(float4) * leaf_size);
    node->leaf.count = n;

    for (int i = 0; i < n; ++i) {
      node->leaf.head_addr[i] = data[low + i];
    }
  } else {
    const int mid = (low + high) / 2;
    const int axis = depth % 4;

    nth_element(data + low, data + mid, data + high, axis);

    node->type = BRANCH;
    node->branch.point = data[mid];
    node->branch.left_child = BuildTree(data, low, mid, leaf_size, depth + 1);
    node->branch.right_child =
        BuildTree(data, mid + 1, high, leaf_size, depth + 1);
  }

  return node;
}

inline Direction flip_dir(const Direction dir) { return 1 - dir; }

inline TreeNode* get_child(const TreeNode* cur, const Direction dir) {
  return dir == LEFT ? cur->branch.left_child : cur->branch.right_child;
}

void traverse_recursive(const TreeNode* cur, const float4 q, float* my_min,
                        const int depth) {
  if (cur->type == LEAF) {
    for (int i = 0; i < depth; ++i) putchar('-');
    printf("%d:%d\n", cur->uid, cur->leaf.count);

    for (int i = 0; i < cur->leaf.count; ++i) {
      const float4 p = cur->leaf.head_addr[i];

      const float dist = kernel_func_4f(p, q);
      *my_min = fminf(*my_min, dist);
    }
  } else {
    for (int i = 0; i < depth; ++i) putchar('-');
    printf("%d\n", cur->uid);

    reduce_element(cur->branch.point, q, my_min);

    const int axis = depth % 4;
    const float train = get_dim(axis, cur->branch.point);
    const float q_value = get_dim(axis, q);
    const Direction dir = q_value < train ? LEFT : RIGHT;

    traverse_recursive(get_child(cur, dir), q, my_min, depth + 1);

    const float diff = kernel_func_1f(q_value, train);
    if (diff < *my_min) {
      traverse_recursive(get_child(cur, flip_dir(dir)), q, my_min, depth + 1);
    }
  }
}

int main(void) {
  const int leaf_size = 32;

  float4 in_data[N];
  for (int i = 0; i < N; ++i) in_data[i] = generate_random_float4();

  const TreeNode* root = BuildTree(in_data, 0, N, leaf_size, 0);
  printf("num_node_used: %d\n", next_node);

  const float4 q = generate_random_float4();
  const float gt = ground_truth(in_data, q);
  printf("gt: %f\n", gt);

  float my_min = FLT_MAX;
  traverse_recursive(root, q, &my_min, 0);
  printf("my_min: %f\n", my_min);

  return 0;
}
