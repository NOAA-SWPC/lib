/* lapack_lls.h
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

#ifndef INCLUDED_lapack_lls_h
#define INCLUDED_lapack_lls_h

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Prototypes
 */

int lls(const gsl_matrix *A, const gsl_vector *c, gsl_vector *x);
int wlls(const gsl_matrix *A, const gsl_vector *w, const gsl_vector *c,
         gsl_vector *x);

#endif /* INCLUDED_lapack_lls_h */
