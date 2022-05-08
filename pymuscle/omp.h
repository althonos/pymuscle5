#pragma once

typedef int omp_lock_t;

static inline int omp_get_thread_num() { return 0; }
static inline int omp_get_max_threads() { return 1; }
static inline void omp_init_lock(omp_lock_t* lock) {}
static inline void omp_set_lock(omp_lock_t* lock) {}
static inline void omp_unset_lock(omp_lock_t* lock) {}
