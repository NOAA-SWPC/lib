/*
 * ema.h
 */

#ifndef INCLUDED_ema_h
#define INCLUDED_ema_h

int ema(const double alpha, double y[], const size_t n);
int ema_reverse(const double alpha, double y[], const size_t n);

#endif /* INCLUDED_ema_h */
