/* lse.h
 * 
 * Copyright (C) 2007 Patrick Alken
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

#ifndef INCLUDED_lse_h
#define INCLUDED_lse_h

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Prototypes
 */

int lse(const gsl_matrix *A, const gsl_matrix *B, const gsl_vector *c,
        const gsl_vector *d, gsl_vector *x);
int wlse(const gsl_matrix *A, const gsl_matrix *B, const gsl_vector *c,
         const gsl_vector *d, const gsl_vector *w, gsl_vector *x);

#endif /* INCLUDED_lse_h */
