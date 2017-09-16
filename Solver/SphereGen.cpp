/***************************************************************************
*                                                                          *
*   Copyright (c) 2017                                                     *
*   FastFieldSolvers S.R.L.  http://www.fastfieldsolvers.com               *
*                                                                          *
*   This program is free software; you can redistribute it and/or modify   *
*   it under the terms of the GNU Lesser General Public License (LGPL)     *
*   as published by the Free Software Foundation; either version 2 of      *
*   the License, or (at your option) any later version.                    *
*   for detail see the LICENCE text file.                                  *
*                                                                          *
*   This program is distributed in the hope that it will be useful,        *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
*   GNU Library General Public License for more details.                   *
*                                                                          *
*   You should have received a copy of the GNU Library General Public      *
*   License along with this program; if not, write to the Free Software    *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307   *
*   USA                                                                    *
*                                                                          *
***************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define X .525731112119133606
#define Z .850650808352039932

static double vdata[12][3] = {
	{-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},
	{0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
	{Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}  };

static int tindices[20][3] = {
	{0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
	{8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
	{7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
	{6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}  };

FILE *fp;
double radius;

void drawtriangle(double *v1, double *v2, double *v3)
{
	fprintf(fp, "T sphere  %f %f %f  %f %f %f  %f %f %f\n", v1[0] * radius, v1[1] * radius, v1[2] * radius,
					v2[0] * radius, v2[1] * radius, v2[2] * radius,
					v3[0] * radius, v3[1] * radius, v3[2] * radius);
}

void normalize(double v[3])
{
	double d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

	if (d == 0.0) {
		fprintf(stderr, "ERROR: Zero length vector!\n");
		return;
	}

	v[0] /= d;
	v[1] /= d;
	v[2] /= d;
}

void subdivide(double *v1, double *v2, double *v3, long depth)
{
	double v12[3], v23[3], v31[3];
	int i;

	if (depth == 0) {
		drawtriangle(v1, v2, v3);
		return;
	}

	for (i = 0; i < 3; i++) {
		v12[i] = v1[i] + v2[i];
		v23[i] = v2[i] + v3[i];
		v31[i] = v3[i] + v1[i];
	}

	normalize(v12);
	normalize(v23);
	normalize(v31);

	subdivide(v1, v12, v31, depth-1);
	subdivide(v2, v23, v12, depth-1);
	subdivide(v3, v31, v23, depth-1);
	subdivide(v12, v23, v31, depth-1);
}

int gensphere(double sradius, long depth, FILE *sfp)
{
//	char *endnumber;
	int i;

	radius = sradius;
	fp = sfp;

	if(fp == NULL)
		return -1;

	fprintf(fp, "0 Sphere with radius %f, refined up to depth %d\n", radius, (int)depth);

	for (i = 0; i < 20; i++) {
		subdivide(&vdata[tindices[i][0]][0],
					&vdata[tindices[i][1]][0],
					&vdata[tindices[i][2]][0], depth);
	}

	return 0;
}

