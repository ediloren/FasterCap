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


// file grouping all basic function tests
// E. Di Lorenzo, 2013/02/01

#include "test.h"

#include <math.h>

#include "FasterCapGlobal.h"

#include "Solver/Potential.h"

void test_pot2D()
{

	CPotential potential;
	double result, res[20], h, enplus, enminus, distmod;
	C2DVector point, normal, dist;
	C2DVector_float vertexes1[2], center;


    LogMsg("Potential 2D\n");

	vertexes1[0].x = -1.0;
	vertexes1[0].y = 0;
	vertexes1[1].x = 1.0;
	vertexes1[1].y = 0.0;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 0.0;
	point.y = 100.0;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1);
	// verify against 1st approx
	res[0] = -log(sqrt((point.x-center.x)*(point.x-center.x)+(point.y-center.y)*(point.y-center.y)));

    LogMsg("Result 1: %g, res[0]: %g\n", result, res[0]);

	// numerical result
	//for(i=1; i<=17; i++) {
	//	res[i-1] = potential.PotentialNumerical(point, vertexes1, i);
	//}

	vertexes1[0].x = -1.0;
	vertexes1[0].y = 0;
	vertexes1[1].x = 1.0;
	vertexes1[1].y = 0.0;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 0.0;
	point.y = 10.0;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1);
	// verify against 1st approx
	res[0] = -log(sqrt((point.x-center.x)*(point.x-center.x)+(point.y-center.y)*(point.y-center.y)));

    LogMsg("Result 2: %g, res[0]: %g\n", result, res[0]);

	vertexes1[0].x = 0.0;
	vertexes1[0].y = -1.0;
	vertexes1[1].x = 0.0;
	vertexes1[1].y = 1.0;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 100.0;
	point.y = 0.0;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1);
	// verify against 1st approx
	res[0] = -log(sqrt((point.x-center.x)*(point.x-center.x)+(point.y-center.y)*(point.y-center.y)));

    LogMsg("Result 3: %g, res[0]: %g\n", result, res[0]);

	vertexes1[0].x = 0.5f;
	vertexes1[0].y = -1.3f;
	vertexes1[1].x = 0.9f;
	vertexes1[1].y = 2.0f;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 28.0;
	point.y = 72.0;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1);
	// verify against 1st approx
	res[0] = -log(sqrt((point.x-center.x)*(point.x-center.x)+(point.y-center.y)*(point.y-center.y)));

    LogMsg("Result 4: %g, res[0]: %g\n", result, res[0]);


	vertexes1[0].x = 0.5f;
	vertexes1[0].y = -1.3f;
	vertexes1[1].x = 0.9f;
	vertexes1[1].y = 2.0f;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = -12.0;
	point.y = -5.0;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1);
	// verify against 1st approx
	res[0] = -log(sqrt((point.x-center.x)*(point.x-center.x)+(point.y-center.y)*(point.y-center.y)));

    LogMsg("Result 5: %g, res[0]: %g\n", result, res[0]);

	vertexes1[0].x = -1.0;
	vertexes1[0].y = 0.0;
	vertexes1[1].x = 1.0;
	vertexes1[1].y = 0.0;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 0.0;
	point.y = 0.1;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1) / TWO_PI_TIMES_E0;
	// verify against 1st approx
	res[0] = -log(sqrt((point.x-center.x)*(point.x-center.x)+(point.y-center.y)*(point.y-center.y))) / TWO_PI_TIMES_E0;

    LogMsg("Actual pot for pplate, mutual: %g, res[0]: %g\n", result, res[0]);

	vertexes1[0].x = -1.0;
	vertexes1[0].y = 0.0;
	vertexes1[1].x = 1.0;
	vertexes1[1].y = 0.0;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 0.0;
	point.y = 0.0;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1) / TWO_PI_TIMES_E0;

    LogMsg("Actual pot for pplate, self: %g\n", result);

	vertexes1[0].x = -1.0;
	vertexes1[0].y = 0.0;
	vertexes1[1].x = 0.0;
	vertexes1[1].y = 0.0;
	point.x = -0.5;
	point.y = 1.0;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1) / TWO_PI_TIMES_E0;

    LogMsg("Mutual cross pot for pplate 1: %g\n", result);

	vertexes1[0].x = -1.0;
	vertexes1[0].y = 0.0;
	vertexes1[1].x = 0.0;
	vertexes1[1].y = 0.0;
	point.x = 0.5;
	point.y = 1.0;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1) / TWO_PI_TIMES_E0;

    LogMsg("Mutual cross pot for pplate 2: %g\n", result);

	vertexes1[0].x = -1.0;
	vertexes1[0].y = 0.0;
	vertexes1[1].x = 0.0;
	vertexes1[1].y = 0.0;
	point.x = 0.5;
	point.y = 0.0;

	// 'result' should be around ???
	result = potential.PotentialOpt(point, vertexes1) / TWO_PI_TIMES_E0;

    LogMsg("Mutual self cross pot for pplate 2: %g\n", result);

    //
    // electric field
    //

    LogMsg("Electric field 2D\n");

	vertexes1[0].x = -1.0;
	vertexes1[0].y = 0;
	vertexes1[1].x = 1.0;
	vertexes1[1].y = 0.0;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 0.0;
	point.y = 100.0;
	normal.x = 0.0;
	normal.y = 1.0;
	h = Mod(vertexes1[1]-vertexes1[0]) / 20.0;

	// 'result' should be around ???
	result = potential.EnField(point, vertexes1, normal);
	// verify against divided differences
	enplus = potential.PotentialOpt(point + normal * h, vertexes1);
	enminus = potential.PotentialOpt(point - normal * h, vertexes1);
	res[0] = -(enplus - enminus) / (2*h);
	// and against 1st approx
	dist = point-center;
	distmod = Mod(dist);
	res[1] = DotProd(dist, normal) / (distmod * distmod);

    LogMsg("Result 1: %g, res[0]: %g, res[1]: %g\n", result, res[0], res[1]);


	vertexes1[0].x = -1.0;
	vertexes1[0].y = 0;
	vertexes1[1].x = 1.0;
	vertexes1[1].y = 0.0;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 0.0;
	point.y = 10.0;
	normal.x = 0.0;
	normal.y = 1.0;
	h = Mod(vertexes1[1]-vertexes1[0]) / 20.0;

	// 'result' should be around ???
	result = potential.EnField(point, vertexes1, normal);
	// verify against divided differences
	enplus = potential.PotentialOpt(point + normal * h, vertexes1);
	enminus = potential.PotentialOpt(point - normal * h, vertexes1);
	res[0] = -(enplus - enminus) / (2*h);
	dist = point-center;
	distmod = Mod(dist);
	res[1] = DotProd(dist, normal) / (distmod * distmod);

    LogMsg("Result 2: %g, res[0]: %g, res[1]: %g\n", result, res[0], res[1]);


	vertexes1[0].x = 0.0;
	vertexes1[0].y = -1.0;
	vertexes1[1].x = 0.0;
	vertexes1[1].y = 1.0;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 100.0;
	point.y = 0.0;
	normal.x = 0.0;
	normal.y = 1.0;
	h = Mod(vertexes1[1]-vertexes1[0]) / 20.0;

	// 'result' should be around ???
	result = potential.EnField(point, vertexes1, normal);
	// verify against divided differences
	enplus = potential.PotentialOpt(point + normal * h, vertexes1);
	enminus = potential.PotentialOpt(point - normal * h, vertexes1);
	res[0] = -(enplus - enminus) / (2*h);
	dist = point-center;
	distmod = Mod(dist);
	res[1] = DotProd(dist, normal) / (distmod * distmod);

    LogMsg("Result 3: %g, res[0]: %g, res[1]: %g\n", result, res[0], res[1]);


	vertexes1[0].x = 0.5f;
	vertexes1[0].y = -1.3f;
	vertexes1[1].x = 0.9f;
	vertexes1[1].y = 2.0f;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = 28.0;
	point.y = 72.0;
	normal.x = 0.0;
	normal.y = 1.0;
	h = Mod(vertexes1[1]-vertexes1[0]) / 20.0;

	// 'result' should be around ???
	result = potential.EnField(point, vertexes1, normal);
	// verify against divided differences
	enplus = potential.PotentialOpt(point + normal * h, vertexes1);
	enminus = potential.PotentialOpt(point - normal * h, vertexes1);
	res[0] = -(enplus - enminus) / (2*h);
	dist = point-center;
	distmod = Mod(dist);
	res[1] = DotProd(dist, normal) / (distmod * distmod);

    LogMsg("Result 4: %g, res[0]: %g, res[1]: %g\n", result, res[0], res[1]);


	vertexes1[0].x = 0.5f;
	vertexes1[0].y = -1.3f;
	vertexes1[1].x = 0.9f;
	vertexes1[1].y = 2.0f;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = -12.0;
	point.y = -5.0;
	normal.x = 0.0;
	normal.y = 1.0;
	h = Mod(vertexes1[1]-vertexes1[0]) / 20.0;

	// 'result' should be around ???
	result = potential.EnField(point, vertexes1, normal);
	// verify against divided differences
	enplus = potential.PotentialOpt(point + normal * h, vertexes1);
	enminus = potential.PotentialOpt(point - normal * h, vertexes1);
	res[0] = -(enplus - enminus) / (2*h);
	dist = point-center;
	distmod = Mod(dist);
	res[1] = DotProd(dist, normal) / (distmod * distmod);

    LogMsg("Result 5: %g, res[0]: %g, res[1]: %g\n", result, res[0], res[1]);


	vertexes1[0].x = 0.5f;
	vertexes1[0].y = -1.3f;
	vertexes1[1].x = 0.9f;
	vertexes1[1].y = 2.0f;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = -12.0;
	point.y = -5.0;
	normal.x = 1.0;
	normal.y = 0.0;
	h = Mod(vertexes1[1]-vertexes1[0]) / 20.0;

	// 'result' should be around ???
	result = potential.EnField(point, vertexes1, normal);
	// verify against divided differences
	enplus = potential.PotentialOpt(point + normal * h, vertexes1);
	enminus = potential.PotentialOpt(point - normal * h, vertexes1);
	res[0] = -(enplus - enminus) / (2*h);
	dist = point-center;
	distmod = Mod(dist);
	res[1] = DotProd(dist, normal) / (distmod * distmod);

    LogMsg("Result 5.1: %g, res[0]: %g, res[1]: %g\n", result, res[0], res[1]);


	vertexes1[0].x = 0.5f;
	vertexes1[0].y = -1.3f;
	vertexes1[1].x = 0.9f;
	vertexes1[1].y = 2.0f;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = -12.0;
	point.y = -5.0;
	normal.x = 1.4;
	normal.y = 1.4;
	normal.Normalize();
	h = Mod(vertexes1[1]-vertexes1[0]) / 20.0;

	// 'result' should be around ???
	result = potential.EnField(point, vertexes1, normal);
	// verify against divided differences
	enplus = potential.PotentialOpt(point + normal * h, vertexes1);
	enminus = potential.PotentialOpt(point - normal * h, vertexes1);
	res[0] = -(enplus - enminus) / (2*h);
	dist = point-center;
	distmod = Mod(dist);
	res[1] = DotProd(dist, normal) / (distmod * distmod);

    LogMsg("Result 5.2: %g, res[0]: %g, res[1]: %g\n", result, res[0], res[1]);


	vertexes1[0].x = 0.5f;
	vertexes1[0].y = -1.3f;
	vertexes1[1].x = 0.9f;
	vertexes1[1].y = 2.0f;
	center = (vertexes1[1] + vertexes1[0]) / 2.0;
	point.x = -12.0;
	point.y = -5.0;
	normal.x = -3.4;
	normal.y = 2.8;
	normal.Normalize();
	h = Mod(vertexes1[1]-vertexes1[0]) / 20.0;

	// 'result' should be around ???
	result = potential.EnField(point, vertexes1, normal);
	// verify against divided differences
	enplus = potential.PotentialOpt(point + normal * h, vertexes1);
	enminus = potential.PotentialOpt(point - normal * h, vertexes1);
	res[0] = -(enplus - enminus) / (2*h);
	dist = point-center;
	distmod = Mod(dist);
	res[1] = DotProd(dist, normal) / (distmod * distmod);

    LogMsg("Result 5.3: %g, res[0]: %g, res[1]: %g\n", result, res[0], res[1]);

}
