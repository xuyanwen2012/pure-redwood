#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

// ---------------------------------------------------------------------------
//  App Constants
// ---------------------------------------------------------------------------

#ifdef REDWOOD_DEBUG
#define DEBUG_PRINT(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...) \
  do {                        \
  } while (0)
#endif

enum {
  N = 10240,
  M = 8,
  MAX_NODES = 2048,
  MAX_STACK_SIZE = 32,

  // how many elements can be processed by the FPGA
  DUET_LEAF_SIZE = 32,
  NUM_EXECUTORS = 2,
};

// ---------------------------------------------------------------------------
//  Utils / Tree Node
// ---------------------------------------------------------------------------

typedef struct float4 {
  float x, y, z, w;
} float4;

float4 generate_random_float4(void) {
  float4 random_float4;
  random_float4.x = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.y = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.z = (float)rand() / RAND_MAX * 1000.0f;
  random_float4.w = (float)rand() / RAND_MAX * 1000.0f;
  return random_float4;
}

typedef enum { LEFT = 0, RIGHT = 1 } Direction;

typedef struct Node {
  int axis;
  int left;
  int right;
  float4 point;
} Node;

typedef struct Range {
  int low;
  int high;
} Range;

// The actual data
float4 in_data[N];
Node nodes[MAX_NODES];
Range ranges[MAX_NODES];

int next_node = 0;

int new_node(void) {
  const int cur = next_node;
  nodes[cur].left = -1;
  nodes[cur].right = -1;
  ranges[cur].low = -1;
  ranges[cur].high = -1;
  ++next_node;
  return cur;
}

bool is_leaf(const int idx) {
  return nodes[idx].left == -1 && nodes[idx].right == -1;
}

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

float kernel_func_4f(const float4 p, const float4 q) {
  const float dx = p.x - q.x;
  const float dy = p.y - q.y;
  const float dz = p.z - q.z;
  const float dw = p.w - q.w;
  return sqrtf(dx * dx + dy * dy + dz * dz + dw * dw);
}

float kernel_func_1f(const float p, const float q) {
  const float dx = p - q;
  return sqrtf(dx * dx);
}

float ground_truth(const float4* data, const float4 q) {
  float gt = FLT_MAX;
  for (int i = 0; i < N; ++i) {
    const float dist = kernel_func_4f(data[i], q);
    gt = fminf(gt, dist);
  }
  return gt;
}

void cpu_emulated_reduce_32(const float4* in_addr, float* out_addr,
                            const float4 q, const int active) {
  float my_min = FLT_MAX;
  for (int i = 0; i < active; ++i) {
    const float dist = kernel_func_4f(in_addr[i], q);
    my_min = fminf(my_min, dist);
  }
  *out_addr = fminf(*out_addr, my_min);
}

// Reduce all values between range
void fpga_kernel(const int low, const int high, float* out_addr,
                 const float4 q) {
  const int n = high - low;
  const float4* addr = &in_data[low];
  const int remainder = n % DUET_LEAF_SIZE;

  int i = 0;
  for (; i < n; i += DUET_LEAF_SIZE) {
    cpu_emulated_reduce_32(addr + i, out_addr, q, DUET_LEAF_SIZE);
  }

  if (remainder > 0) {
    cpu_emulated_reduce_32(addr + i, out_addr, q, remainder);
  }
}

// ---------------------------------------------------------------------------
//  C++ STD Algorithm (Modified)
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
    // divide and conquer, ordering partition containing Nth
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
//  KDT Tree Data Structure
// ---------------------------------------------------------------------------

int BuildTree(float4* data, const int low, const int high, const int leaf_size,
              const int depth) {
  const int cur = new_node();

  const int n = high - low;
  if (n <= leaf_size) {
    ranges[cur].low = low;
    ranges[cur].high = high;
  } else {
    const int mid = (low + high) / 2;
    const int axis = depth % 4;

    nth_element(data + low, data + mid, data + high, axis);

    nodes[cur].axis = axis;
    nodes[cur].point = data[mid];
    nodes[cur].left = BuildTree(data, low, mid, leaf_size, depth + 1);
    nodes[cur].right = BuildTree(data, mid + 1, high, leaf_size, depth + 1);

    ranges[cur].low = ranges[nodes[cur].left].low;
    ranges[cur].high = ranges[nodes[cur].right].high;
  }

  return cur;
}

Direction flip_dir(const Direction dir) { return 1 - dir; }

int get_child(const int cur, const Direction dir) {
  return dir == LEFT ? nodes[cur].left : nodes[cur].right;
}

void dfs(const int cur, const int depth) {
  if (is_leaf(cur)) {
    for (int i = 0; i < depth; ++i) putchar('-');
    DEBUG_PRINT("[%d]\t[%d, %d)\n", cur, ranges[cur].low, ranges[cur].high);
  } else {
    for (int i = 0; i < depth; ++i) putchar('-');
    DEBUG_PRINT("[%d]\tleft: %d\tright: %d\t[%d, %d)\n", cur, nodes[cur].left,
                nodes[cur].right, ranges[cur].low, ranges[cur].high);

    dfs(get_child(cur, LEFT), depth + 1);
    dfs(get_child(cur, RIGHT), depth + 1);
  }
}

void traverse_recursive(const int cur, const float4 q, float* my_min,
                        const int depth) {
  if (is_leaf(cur)) {
    for (int i = 0; i < depth; ++i) putchar('-');
    DEBUG_PRINT("[%d]\t[%d, %d)\n", cur, ranges[cur].low, ranges[cur].high);

    // Reduce at leaf node, accelerated by FPGA
    fpga_kernel(ranges[cur].low, ranges[cur].high, my_min, q);
  } else {
    for (int i = 0; i < depth; ++i) putchar('-');
    DEBUG_PRINT("[%d]\tleft: %d\tright: %d\t[%d, %d)\n", cur, nodes[cur].left,
                nodes[cur].right, ranges[cur].low, ranges[cur].high);

    // Reduce at branch node
    const float dist = kernel_func_4f(nodes[cur].point, q);
    *my_min = fminf(*my_min, dist);

    const int axis = depth % 4;
    const float train = get_dim(axis, nodes[cur].point);
    const float q_value = get_dim(axis, q);
    const Direction dir = q_value < train ? LEFT : RIGHT;

    traverse_recursive(get_child(cur, dir), q, my_min, depth + 1);

    const float diff = kernel_func_1f(q_value, train);
    if (diff < *my_min) {
      traverse_recursive(get_child(cur, flip_dir(dir)), q, my_min, depth + 1);
    }
  }
}

// ---------------------------------------------------------------------------
//  Redwood Traverser
// ---------------------------------------------------------------------------

// Node traversal stack
typedef struct Fields {
  int cur;
  Direction dir;
  float train;
  float q_value;
} Fields;

typedef enum { Finished, Working } ExecuteStatus;

typedef struct Executor {
  ExecuteStatus state;
  float4 q;
  float my_min;
  int cur_node_stack;
  Fields stack[MAX_STACK_SIZE];
} Executor;

void reset_exe(Executor* exe, const float4 q) {
  exe->state = Finished;
  exe->q = q;
  exe->my_min = FLT_MAX;    // result_set->Reset();
  exe->cur_node_stack = 0;  // stack_.clear();
}

void init_exe(Executor* exe) {
  for (int j = 0; j < MAX_STACK_SIZE; ++j) {
    exe->stack[j].cur = -1;
    exe->stack[j].dir = LEFT;
    exe->stack[j].train = 0.0f;
    exe->stack[j].q_value = 0.0f;
  }
}

int push_node_stack(Executor* exe, const int cur, const Direction dir,
                    const float train, const float q_value) {
  ++exe->cur_node_stack;
  if (exe->cur_node_stack < MAX_STACK_SIZE) {
    exe->stack[exe->cur_node_stack].cur = cur;
    exe->stack[exe->cur_node_stack].dir = dir;
    exe->stack[exe->cur_node_stack].train = train;
    exe->stack[exe->cur_node_stack].q_value = q_value;
    return 0;
  }
  printf("executor stack overflow!!\n");
  exit(EXIT_FAILURE);
}

Fields pop_node_stack(Executor* exe) {
  return exe->stack[exe->cur_node_stack--];
}

bool node_stack_empty(const Executor* exe) { return exe->cur_node_stack == 0; }

void execute(Executor* exe) {
  if (exe->state == Working) {
    goto my_resume_point;
  }

  exe->state = Working;
  exe->cur_node_stack = 0;

  int cur = 0;

  while (cur != -1 || !node_stack_empty(exe)) {
    while (cur != -1) {
      if (is_leaf(cur)) {
        DEBUG_PRINT("(%p) [%d]\t[%d, %d)\n", (void*)exe, cur, ranges[cur].low,
                    ranges[cur].high);

        fpga_kernel(ranges[cur].low, ranges[cur].high, &exe->my_min, exe->q);

        // **** Coroutine Reuturn (API) ****
        return;
      my_resume_point:
        // ****************************

        cur = -1;
        continue;
      }

      DEBUG_PRINT("(%p) [%d]\tleft: %d\tright: %d\t[%d, %d)\n", (void*)exe, cur,
                  nodes[cur].left, nodes[cur].right, ranges[cur].low,
                  ranges[cur].high);

      // Reduce at branch node
      const float dist = kernel_func_4f(nodes[cur].point, exe->q);
      exe->my_min = fminf(exe->my_min, dist);

      const int axis = nodes[cur].axis;
      const float train = get_dim(axis, nodes[cur].point);
      const float q_value = get_dim(axis, exe->q);
      const Direction dir = q_value < train ? LEFT : RIGHT;

      // Recursion 1
      push_node_stack(exe, cur, dir, train, q_value);
      cur = get_child(cur, dir);
    }

    if (!node_stack_empty(exe)) {
      // pop
      const Fields last = pop_node_stack(exe);

      const float diff = kernel_func_1f(last.q_value, last.train);
      if (diff < exe->my_min) {
        // Recursion 2
        cur = get_child(last.cur, flip_dir(last.dir));
      }
    }
  }
  // Done traversals
  exe->state = Finished;
}

void start_new_execute(Executor* exe, const float4 q) {
  reset_exe(exe, q);
  execute(exe);
}

int main(void) {
  const int leaf_size = 32;
  assert(leaf_size >= DUET_LEAF_SIZE);

  for (int i = 0; i < N; ++i) in_data[i] = generate_random_float4();

  BuildTree(in_data, 0, N, leaf_size, 0);

  DEBUG_PRINT("num_node_used: %d\n", next_node);

  Executor exes[NUM_EXECUTORS];
  for (int i = 0; i < NUM_EXECUTORS; ++i) {
    init_exe(&exes[i]);
  }

  int processed = 0;
  while (processed < M) {
    for (int i = 0; i < NUM_EXECUTORS;) {
      if (exes[i].state == Finished) {
        // pick a new task
        const float4 q = generate_random_float4();

        ++processed;
        start_new_execute(&exes[i], q);
        ++i;
      } else {
        execute(&exes[i]);
        if (exes[i].state == Finished) {
          printf("result: %f\n", exes[i].my_min);
        } else {
          ++i;
        }
      }
    }
  }

  return EXIT_SUCCESS;
}
