// minimal test harness
#pragma once
#include <cstdio>
#include <cstdlib>

static int g_fail = 0, g_checks = 0;
#define CHECK(cond)                                                              \
  do {                                                                           \
    ++g_checks;                                                                  \
    if (!(cond)) {                                                               \
      ++g_fail;                                                                  \
      if (g_fail <= 20) fprintf(stderr, "FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); \
    }                                                                            \
  } while (0)
#define CHECKM(cond, ...)                                                        \
  do {                                                                           \
    ++g_checks;                                                                  \
    if (!(cond)) {                                                               \
      ++g_fail;                                                                  \
      if (g_fail <= 20) { fprintf(stderr, "FAIL %s:%d  %s  ", __FILE__, __LINE__, #cond); fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); } \
    }                                                                            \
  } while (0)
static int finish(const char *name) {
  printf("%s: %d checks, %d failures\n", name, g_checks, g_fail);
  return g_fail ? 1 : 0;
}
