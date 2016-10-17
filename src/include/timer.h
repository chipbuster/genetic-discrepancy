#ifndef TIMER_H
#define TIMER_H
#include<sys/time.h>

#define INIT_TIMER(t) struct timeval start_##t, end_##t;
#define START_TIMER(t) gettimeofday(&start_##t, 0);
#define STOP_TIMER(t)  gettimeofday(&end_##t, 0);
#define ELAPSED_TIME(t)  (1000000.0*(end_##t.tv_sec-start_##t.tv_sec) + end_##t.tv_usec-start_##t.tv_usec)/1000000.0
#define PRINT_TIMER(t,s) fprintf(stderr, #t " - Elapsed time for %s: %3.4fs\n", s, ELAPSED_TIME(t));

#endif  // TIMER_H
