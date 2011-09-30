#ifdef TIMING

#include <sys/time.h>

extern struct timezone timing_tz;
extern struct timeval  timing_tv1;
extern struct timeval  timing_tv2;


#define INIT_TIME gettimeofday(&timing_tv1, &timing_tz);

#define PRINT_TIME(a) gettimeofday(&timing_tv2, &timing_tz); \
  fprintf(stderr, "Elapsed time (%s): %d us\n", (a), \
    (int)(((double)(timing_tv2.tv_sec*1e6+timing_tv2.tv_usec) - \
          (double)(timing_tv1.tv_sec*1e6+timing_tv1.tv_usec))));

#else

#define INIT_TIME 
#define PRINT_TIME(a)

#endif
