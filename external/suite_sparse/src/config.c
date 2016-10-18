/* ========================================================================== */
/* === SuiteSparse_config =================================================== */
/* ========================================================================== */

/* SuiteSparse configuration : memory manager and printf functions. */

/* Copyright (c) 2013, Timothy A. Davis.  No licensing restrictions
 * apply to this file or to the SuiteSparse_config directory.
 * Author: Timothy A. Davis.
 */

#include <math.h>
#include <stdlib.h>

#ifndef NPRINT
#include <stdio.h>
#endif

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif

#ifndef NULL
#define NULL ((void *) 0)
#endif

#include "ElSuiteSparse/config.h"

/* -------------------------------------------------------------------------- */
/* SuiteSparse_config : a global extern struct */
/* -------------------------------------------------------------------------- */

/* The SuiteSparse_config struct is available to all SuiteSparse functions and
    to all applications that use those functions.  It must be modified with
    care, particularly in a multithreaded context.  Normally, the application
    will initialize this object once, via SuiteSparse_start, possibily followed
    by application-specific modifications if the applications wants to use
    alternative memory manager functions.

    The user can redefine these global pointers at run-time to change the
    memory manager and printf function used by SuiteSparse.

    If -DNMALLOC is defined at compile-time, then no memory-manager is
    specified.  You must define them at run-time, after calling
    SuiteSparse_start.

    If -DPRINT is defined a compile time, then printf is disabled, and
    SuiteSparse will not use printf.
 */

struct ElSuiteSparse_config_struct ElSuiteSparse_config =
{

    /* memory management functions */
    #ifndef NMALLOC
        #ifdef MATLAB_MEX_FILE
            /* MATLAB mexFunction: */
            mxMalloc, mxCalloc, mxRealloc, mxFree,
        #else
            /* standard ANSI C: */
            malloc, calloc, realloc, free,
        #endif
    #else
        /* no memory manager defined; you must define one at run-time: */
        NULL, NULL, NULL, NULL,
    #endif

    /* printf function */
    #ifndef NPRINT
        #ifdef MATLAB_MEX_FILE
            /* MATLAB mexFunction: */
            mexPrintf,
        #else
            /* standard ANSI C: */
            printf,
        #endif
    #else
        /* printf is disabled */
        NULL,
    #endif
} ;

/* -------------------------------------------------------------------------- */
/* SuiteSparse_start */
/* -------------------------------------------------------------------------- */

/* All applications that use SuiteSparse should call SuiteSparse_start prior
   to using any SuiteSparse function.  Only a single thread should call this
   function, in a multithreaded application.  Currently, this function is
   optional, since all this function currently does is to set the four memory
   function pointers to NULL (which tells SuiteSparse to use the default
   functions).  In a multi- threaded application, only a single thread should
   call this function.

   Future releases of SuiteSparse might enforce a requirement that
   SuiteSparse_start be called prior to calling any SuiteSparse function.
 */

void ElSuiteSparse_start ( void )
{

    /* memory management functions */
    #ifndef NMALLOC
        #ifdef MATLAB_MEX_FILE
            /* MATLAB mexFunction: */
            ElSuiteSparse_config.malloc_func  = mxMalloc ;
            ElSuiteSparse_config.calloc_func  = mxCalloc ;
            ElSuiteSparse_config.realloc_func = mxRealloc ;
            ElSuiteSparse_config.free_func    = mxFree ;
        #else
            /* standard ANSI C: */
            ElSuiteSparse_config.malloc_func  = malloc ;
            ElSuiteSparse_config.calloc_func  = calloc ;
            ElSuiteSparse_config.realloc_func = realloc ;
            ElSuiteSparse_config.free_func    = free ;
        #endif
    #else
        /* no memory manager defined; you must define one after calling
           SuiteSparse_start */
        ElSuiteSparse_config.malloc_func  = NULL ;
        ElSuiteSparse_config.calloc_func  = NULL ;
        ElSuiteSparse_config.realloc_func = NULL ;
        ElSuiteSparse_config.free_func    = NULL ;
    #endif

    /* printf function */
    #ifndef NPRINT
        #ifdef MATLAB_MEX_FILE
            /* MATLAB mexFunction: */
            ElSuiteSparse_config.printf_func = mexPrintf ;
        #else
            /* standard ANSI C: */
            ElSuiteSparse_config.printf_func = printf ;
        #endif
    #else
        /* printf is disabled */
        ElSuiteSparse_config.printf_func = NULL ;
    #endif
}

/* -------------------------------------------------------------------------- */
/* SuiteSparse_finish */
/* -------------------------------------------------------------------------- */

/* This currently does nothing, but in the future, applications should call
   SuiteSparse_start before calling any SuiteSparse function, and then
   SuiteSparse_finish after calling the last SuiteSparse function, just before
   exiting.  In a multithreaded application, only a single thread should call
   this function.

   Future releases of SuiteSparse might use this function for any
   SuiteSparse-wide cleanup operations or finalization of statistics.
 */

void ElSuiteSparse_finish ( void )
{
    /* do nothing */ ;
}

/* -------------------------------------------------------------------------- */
/* SuiteSparse_malloc: malloc wrapper */
/* -------------------------------------------------------------------------- */

void *ElSuiteSparse_malloc    /* pointer to allocated block of memory */
(
    size_t nitems,          /* number of items to malloc */
    size_t size_of_item     /* sizeof each item */
)
{
    void *p ;
    size_t size ;
    if (nitems < 1) nitems = 1 ;
    if (size_of_item < 1) size_of_item = 1 ;
    size = nitems * size_of_item  ;

    if (size != ((double) nitems) * size_of_item)
    {
        /* size_t overflow */
        p = NULL ;
    }
    else
    {
        p = (void *) (ElSuiteSparse_config.malloc_func) (size) ;
    }
    return (p) ;
}


/* -------------------------------------------------------------------------- */
/* SuiteSparse_calloc: calloc wrapper */
/* -------------------------------------------------------------------------- */

void *ElSuiteSparse_calloc    /* pointer to allocated block of memory */
(
    size_t nitems,          /* number of items to calloc */
    size_t size_of_item     /* sizeof each item */
)
{
    void *p ;
    size_t size ;
    if (nitems < 1) nitems = 1 ;
    if (size_of_item < 1) size_of_item = 1 ;
    size = nitems * size_of_item  ;

    if (size != ((double) nitems) * size_of_item)
    {
        /* size_t overflow */
        p = NULL ;
    }
    else
    {
        p = (void *) (ElSuiteSparse_config.calloc_func) (nitems, size_of_item) ;
    }
    return (p) ;
}

/* -------------------------------------------------------------------------- */
/* SuiteSparse_realloc: realloc wrapper */
/* -------------------------------------------------------------------------- */

/* If p is non-NULL on input, it points to a previously allocated object of
   size nitems_old * size_of_item.  The object is reallocated to be of size
   nitems_new * size_of_item.  If p is NULL on input, then a new object of that
   size is allocated.  On success, a pointer to the new object is returned,
   and ok is returned as 1.  If the allocation fails, ok is set to 0 and a
   pointer to the old (unmodified) object is returned.
 */

void *ElSuiteSparse_realloc   /* pointer to reallocated block of memory, or
                               to original block if the realloc failed. */
(
    size_t nitems_new,      /* new number of items in the object */
    size_t nitems_old,      /* old number of items in the object */
    size_t size_of_item,    /* sizeof each item */
    void *p,                /* old object to reallocate */
    int *ok                 /* 1 if successful, 0 otherwise */
)
{
    size_t size ;
    if (nitems_old < 1) nitems_old = 1 ;
    if (nitems_new < 1) nitems_new = 1 ;
    if (size_of_item < 1) size_of_item = 1 ;
    size = nitems_new * size_of_item  ;

    if (size != ((double) nitems_new) * size_of_item)
    {
        /* size_t overflow */
        (*ok) = 0 ;
    }
    else if (p == NULL)
    {
        /* a fresh object is being allocated */
        p = ElSuiteSparse_malloc (nitems_new, size_of_item) ;
        (*ok) = (p != NULL) ;
    }
    else if (nitems_old == nitems_new)
    {
        /* the object does not change; do nothing */
        (*ok) = 1 ;
    }
    else
    {
        /* change the size of the object from nitems_old to nitems_new */
        void *pnew ;
        pnew = (void *) (ElSuiteSparse_config.realloc_func) (p, size) ;
        if (pnew == NULL)
        {
            if (nitems_new < nitems_old)
            {
                /* the attempt to reduce the size of the block failed, but
                   the old block is unchanged.  So pretend to succeed. */
                (*ok) = 1 ;
            }
            else
            {
                /* out of memory */
                (*ok) = 0 ;
            }
        }
        else
        {
            /* success */
            p = pnew ;
            (*ok) = 1 ;
        }
    }
    return (p) ;
}

/* -------------------------------------------------------------------------- */
/* SuiteSparse_free: free wrapper */
/* -------------------------------------------------------------------------- */

void *ElSuiteSparse_free      /* always returns NULL */
(
    void *p                 /* block to free */
)
{
    if (p)
    {
        (ElSuiteSparse_config.free_func) (p) ;
    }
    return (NULL) ;
}


/* -------------------------------------------------------------------------- */
/* SuiteSparse_tic: return current wall clock time */
/* -------------------------------------------------------------------------- */

/* Returns the number of seconds (tic [0]) and nanoseconds (tic [1]) since some
 * unspecified but fixed time in the past.  If no timer is installed, zero is
 * returned.  A scalar double precision value for 'tic' could be used, but this
 * might cause loss of precision because clock_getttime returns the time from
 * some distant time in the past.  Thus, an array of size 2 is used.
 *
 * The timer is enabled by default.  To disable the timer, compile with
 * -DNTIMER.  If enabled on a POSIX C 1993 system, the timer requires linking
 * with the -lrt library.
 *
 * example:
 *
 *      double tic [2], r, s, t ;
 *      SuiteSparse_tic (tic) ;     // start the timer
 *      // do some work A
 *      t = SuiteSparse_toc (tic) ; // t is time for work A, in seconds
 *      // do some work B
 *      s = SuiteSparse_toc (tic) ; // s is time for work A and B, in seconds
 *      SuiteSparse_tic (tic) ;     // restart the timer
 *      // do some work C
 *      r = SuiteSparse_toc (tic) ; // s is time for work C, in seconds
 *
 * A double array of size 2 is used so that this routine can be more easily
 * ported to non-POSIX systems.  The caller does not rely on the POSIX
 * <time.h> include file.
 */

#ifdef EL_SUITESPARSE_TIMER_ENABLED

#include <time.h>

void ElSuiteSparse_tic
(
    double tic [2]      /* output, contents undefined on input */
)
{
    /* POSIX C 1993 timer, requires -librt */
    struct timespec t ;
    clock_gettime (CLOCK_MONOTONIC, &t) ;
    tic [0] = (double) (t.tv_sec) ;
    tic [1] = (double) (t.tv_nsec) ;
}

#else

void ElSuiteSparse_tic
(
    double tic [2]      /* output, contents undefined on input */
)
{
    /* no timer installed */
    tic [0] = 0 ;
    tic [1] = 0 ;
}

#endif


/* -------------------------------------------------------------------------- */
/* SuiteSparse_toc: return time since last tic */
/* -------------------------------------------------------------------------- */

/* Assuming SuiteSparse_tic is accurate to the nanosecond, this function is
 * accurate down to the nanosecond for 2^53 nanoseconds since the last call to
 * SuiteSparse_tic, which is sufficient for SuiteSparse (about 104 days).  If
 * additional accuracy is required, the caller can use two calls to
 * SuiteSparse_tic and do the calculations differently.
 */

double ElSuiteSparse_toc  /* returns time in seconds since last tic */
(
    double tic [2]  /* input, not modified from last call to SuiteSparse_tic */
)
{
    double toc [2] ;
    ElSuiteSparse_tic (toc) ;
    return ((toc [0] - tic [0]) + 1e-9 * (toc [1] - tic [1])) ;
}


/* -------------------------------------------------------------------------- */
/* SuiteSparse_time: return current wallclock time in seconds */
/* -------------------------------------------------------------------------- */

/* This function might not be accurate down to the nanosecond. */

double ElSuiteSparse_time  /* returns current wall clock time in seconds */
(
    void
)
{
    double toc [2] ;
    ElSuiteSparse_tic (toc) ;
    return (toc [0] + 1e-9 * toc [1]) ;
}


/* -------------------------------------------------------------------------- */
/* SuiteSparse_version: return the current version of SuiteSparse */
/* -------------------------------------------------------------------------- */

int ElSuiteSparse_version
(
    int version [3]
)
{
    if (version != NULL)
    {
        version [0] = EL_SUITESPARSE_MAIN_VERSION ;
        version [1] = EL_SUITESPARSE_SUB_VERSION ;
        version [2] = EL_SUITESPARSE_SUBSUB_VERSION ;
    }
    return (EL_SUITESPARSE_VERSION) ;
}
