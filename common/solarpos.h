/* solarpos.h
 * 
 * Copyright (C) 2006, 2007 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_solarpos_h
#define INCLUDED_solarpos_h

#include "solpos00.h"

typedef struct
{
  struct posdata pd;
} solarpos_workspace;

/*
 * Prototypes
 */

solarpos_workspace *solarpos_alloc(void);
void solarpos_free(solarpos_workspace *w);
void solarpos_init(solarpos_workspace *w);
int solarpos_calc_zenith(time_t ts, double latitude, double longitude,
                         double *result, solarpos_workspace *w);
int solarpos_calc_sunrs(time_t t, double latitude, double longitude,
                        double *sunrise, double *sunset,
                        solarpos_workspace *w);

#endif /* INCLUDED_solarpos_h */
