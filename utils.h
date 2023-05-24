#pragma once

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