/*
    OpenUniverse 1.0
    Copyright (C) 2000  Raul Alonso <amil@las.es>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/* Macros used to scale down distances/radii so we can reduce the
  jerkiness effect in outter bodies */
#define AU 149597.870

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define RADIUSSCALE(x) ((x)*6.378)
#define DISTCORRECTION(x) ((x) * AU)

/* Need no comment ;-) */
#define DEG2RAD(x) ((x)*M_PI/180.0)
#define RAD2DEG(x) ((x)*180.0/M_PI)
#define DISTANCE(x,y,z) sqrt((x)*(x)+(y)*(y)+(z)*(z))
#define CLAMP( X, MIN, MAX )  ( (X)<(MIN) ? (MIN) : ((X)>(MAX) ? (MAX) : (X)) )
#define CLAMP_SELF(x, mn, mx)  \
   ( (x)<(mn) ? ((x) = (mn)) : ((x)>(mx) ? ((x)=(mx)) : (x)) )

/* Vector macros */

#define INITVECTOR( A, x, y, z)                 \
do {						\
   (A)[0] = (x);                                \
   (A)[1] = (y);                                \
   (A)[2] = (z);                                \
} while (0)


#define COPYVECTOR( A, B )                  \
do {						\
   (A)[0] = (B)[0];                         \
   (A)[1] = (B)[1];                         \
   (A)[2] = (B)[2];                         \
} while (0)


#define SAMEVECTOR( A, B )                  \
( ((A)[0] == (B)[0]) ? (( ((A)[1] == (B)[1]) ? (( ((A)[2] == (B)[2]) ? (1) : 0)) : 0 )) : 0 )
				


#define MODULE(V) sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2])

#define NORMALIZE(V)                            \
do {                                            \
   double d=MODULE(V);                          \
   (V)[0]/=d;                                   \
   (V)[1]/=d;                                   \
   (V)[2]/=d;                                   \
} while (0)


#define ADDVECTORS( A, B, C )                   \
do {						\
      (A)[0] = (B)[0] + (C)[0];                 \
      (A)[1] = (B)[1] + (C)[1];                 \
      (A)[2] = (B)[2] + (C)[2];                 \
} while (0)


#define SUBVECTORS( A, B, C )                   \
do {						\
      (A)[0] = (B)[0] - (C)[0];                 \
      (A)[1] = (B)[1] - (C)[1];                 \
      (A)[2] = (B)[2] - (C)[2];                 \
} while (0)


#define MULTVECTOR( A, B, K )                   \
do {						\
      (A)[0] = (B)[0] * (K);                      \
      (A)[1] = (B)[1] * (K);                      \
      (A)[2] = (B)[2] * (K);                      \
} while (0)


#define DIVVECTOR( A, B, K )                    \
do {                                            \
      (A)[0] = (B)[0] / (K);                      \
      (A)[1] = (B)[1] / (K);                      \
      (A)[2] = (B)[2] / (K);                      \
} while (0)


#define DOTPRODUCT( A, B )  ( (A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2] )

#define CROSSPRODUCT(A, B, C)                         \
do {						\
   double x,y,z;                                \
   x = (B)[1]*(C)[2] - (B)[2]*(C)[1];           \
   y = (B)[2]*(C)[0] - (B)[0]*(C)[2];           \
   z = (B)[0]*(C)[1] - (B)[1]*(C)[0];           \
   (A)[0] = x;                                  \
   (A)[1] = y;                                  \
   (A)[2] = z;                                  \
} while (0)

enum{X,Y,Z};
