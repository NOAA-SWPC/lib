/*
 * MATLAB Compiler: 6.1 (R2015b)
 * Date: Mon Sep 28 15:24:11 2015
 * Arguments: "-B" "macro_default" "-B" "csharedlib:libchaosext" "-W"
 * "lib:libchaosext" "-T" "link:lib" "-d" "libchaosext" "synth_chaos_ext.m" 
 */

#ifndef __libchaosext_h
#define __libchaosext_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libchaosext
#define PUBLIC_libchaosext_C_API __global
#else
#define PUBLIC_libchaosext_C_API /* No import statement needed. */
#endif

#define LIB_libchaosext_C_API PUBLIC_libchaosext_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libchaosext
#define PUBLIC_libchaosext_C_API __declspec(dllexport)
#else
#define PUBLIC_libchaosext_C_API __declspec(dllimport)
#endif

#define LIB_libchaosext_C_API PUBLIC_libchaosext_C_API


#else

#define LIB_libchaosext_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libchaosext_C_API 
#define LIB_libchaosext_C_API /* No special import/export declaration */
#endif

extern LIB_libchaosext_C_API 
bool MW_CALL_CONV libchaosextInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_libchaosext_C_API 
bool MW_CALL_CONV libchaosextInitialize(void);

extern LIB_libchaosext_C_API 
void MW_CALL_CONV libchaosextTerminate(void);



extern LIB_libchaosext_C_API 
void MW_CALL_CONV libchaosextPrintStackTrace(void);

extern LIB_libchaosext_C_API 
bool MW_CALL_CONV mlxSynth_chaos_ext(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                     *prhs[]);



extern LIB_libchaosext_C_API bool MW_CALL_CONV mlfSynth_chaos_ext(int nargout, mxArray** B_ext, mxArray* t, mxArray* r, mxArray* theta, mxArray* phi, mxArray* RC_e, mxArray* RC_i);

#ifdef __cplusplus
}
#endif
#endif
