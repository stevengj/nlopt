/* Copyright (c) 2007-2010 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#include "nlopt-util.h"

/* Simple replacement for the BSD qsort_r function (re-entrant sorting),
   if it is not available.

   (glibc 2.8 included a qsort_r function as well, but totally
   *%&$#-ed things up by gratuitously changing the argument order, in
   such a way as to allow code using the BSD ordering to compile but
   die a flaming death at runtime.  Damn them all to Hell, I'll just
   use my own implementation.)

   (Actually, with glibc 2.3.6 on my Intel Core Duo, my implementation
   below seems to be significantly faster than qsort.  Go figure.)
*/

#ifndef HAVE_QSORT_R_damn_it_use_my_own
/* swap size bytes between a_ and b_ */
static void swap(void *a_, void *b_, size_t size)
{
     if (a_ == b_) return;
     {
          size_t i, nlong = size / sizeof(long);
          long *a = (long *) a_, *b = (long *) b_;
          for (i = 0; i < nlong; ++i) {
               long c = a[i];
               a[i] = b[i];
               b[i] = c;
          }
	  a_ = (void*) (a + nlong);
	  b_ = (void*) (b + nlong);
     }
     {
          size_t i;
          char *a = (char *) a_, *b = (char *) b_;
          size = size % sizeof(long);
          for (i = 0; i < size; ++i) {
               char c = a[i];
               a[i] = b[i];
               b[i] = c;
          }
     }
}
#endif /* HAVE_QSORT_R */

void nlopt_qsort_r(void *base_, size_t nmemb, size_t size, void *thunk,
		   int (*compar)(void *, const void *, const void *))
{
#ifdef HAVE_QSORT_R_damn_it_use_my_own
     /* Even if we could detect glibc vs. BSD by appropriate
	macrology, there is no way to make the calls compatible
	without writing a wrapper for the compar function...screw
	this. */
     qsort_r(base_, nmemb, size, thunk, compar);
#else
     char *base = (char *) base_;
     if (nmemb < 10) { /* use O(nmemb^2) algorithm for small enough nmemb */
	  size_t i, j;
	  for (i = 0; i+1 < nmemb; ++i)
	       for (j = i+1; j < nmemb; ++j)
		    if (compar(thunk, base+i*size, base+j*size) > 0)
			 swap(base+i*size, base+j*size, size);
     }
     else {
	  size_t i, pivot, npart;
	  /* pick median of first/middle/last elements as pivot */
	  {
	       const char *a = base, *b = base + (nmemb/2)*size, 
		    *c = base + (nmemb-1)*size;
	       pivot = compar(thunk,a,b) < 0
		    ? (compar(thunk,b,c) < 0 ? nmemb/2 :
		       (compar(thunk,a,c) < 0 ? nmemb-1 : 0))
		    : (compar(thunk,a,c) < 0 ? 0 :
		       (compar(thunk,b,c) < 0 ? nmemb-1 : nmemb/2));
	  }
	  /* partition array */
	  swap(base + pivot*size, base + (nmemb-1) * size, size);
	  pivot = (nmemb - 1) * size;
	  for (i = npart = 0; i < nmemb-1; ++i)
	       if (compar(thunk, base+i*size, base+pivot) <= 0)
		    swap(base+i*size, base+(npart++)*size, size);
	  swap(base+npart*size, base+pivot, size);
	  /* recursive sort of two partitions */
	  nlopt_qsort_r(base, npart, size, thunk, compar);
	  npart++; /* don't need to sort pivot */
	  nlopt_qsort_r(base+npart*size, nmemb-npart, size, thunk, compar);
     }
#endif /* !HAVE_QSORT_R */
}
