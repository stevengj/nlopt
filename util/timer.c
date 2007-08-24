#include "nlopt-util.h"
#include "config.h"

#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

#if defined(_WIN32) || defined(__WIN32__)
#  include <windows.h>    
#endif

/* return time in seconds since some arbitrary point in the past */
double nlopt_seconds(void)
{
     static int start_inited = 0; /* whether start time has been initialized */
#if defined(HAVE_GETTIMEOFDAY)
     static struct timeval start;
     struct timeval tv;
     if (!start_inited) {
	  start_inited = 1;
	  gettimeofday(&start, NULL);
     }
     gettimeofday(&tv, NULL);
     return (tv.tv_sec - start.tv_sec) + 1.e-6 * (tv.tv_usec - start.tv_sec);
#elif defined(HAVE_TIME)
     return time(NULL);
#elif defined(_WIN32) || defined(__WIN32__)
     static ULONGLONG start;
     FILETIME ft;
     if (!start_inited) {
	  start_inited = 1;
	  GetSystemTimeAsFileTime(&ft);
	  start = (((ULONGLONG) ft.dwHighDateTime) << 32) + ft.dwLowDateTime;
     }
     GetSystemTimeAsFileTime(&ft);
     return 100e-9 * (((((ULONGLONG) ft.dwHighDateTime) << 32) + ft.dwLowDateTime) - start);
#else
     /* use clock() as a fallback... this is somewhat annoying
	because clock() may wrap around with a fairly short period */
     static clock_t start;
     if (!start_inited) {
	  start_inited = 1;
	  start = clock();
     }
     return (clock() - start) * 1.0 / CLOCKS_PER_SEC;
#endif
}

/* number based on time for use as random seed */
unsigned long nlopt_time_seed(void)
{
#if defined(HAVE_GETTIMEOFDAY)
     struct timeval tv;
     gettimeofday(&tv, NULL);
     return (tv.tv_sec ^ tv.tv_usec);
#elif defined(HAVE_TIME)
     return time(NULL);
#elif defined(_WIN32) || defined(__WIN32__)
     FILETIME ft;
     GetSystemTimeAsFileTime(&ft);
     return ft.dwHighDateTime ^ ft.dwLowDateTime;
#else
     return clock();
#endif
}
