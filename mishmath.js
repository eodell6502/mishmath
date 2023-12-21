/*
Copyright 2019 Erich Waidthaler, original authors, and subsequent contributors

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// Note: It would be ideal if I could use strict mode for the entire module,
// and hopefully we will eventually get there, but in the meantime, it's used
// on a function-by-function basis.

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Util functions which need to move to a separate file... TODO

function uniq(arr) {
    if(arr.length < 2)
        return arr;

    const sorted = arr.sort();

    let result = [arr[0]];
    for(var i = 1; i < sorted.length; i++) {
        if(sorted[i - 1].join('') !== sorted[i].join('')) {
            result.push(sorted[i]);
        }
    }
    return result;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


var mishmath = { };


//==============================================================================
// Implements the fast Fisher-Yates shuffle. If the optional inplace argument
// is true, which is the default, arr is shuffled in place. The optional rng
// argument is a random number generation function. If not specified,
// Math.random is used. (Stupid JS trick: rng can be specified without inplace.)
// The shuffled array is returned.
//==============================================================================

mishmath.fisherYatesShuffle = function(arr, inplace, rng) {
    "use strict";

    if(inplace === undefined) {
        inplace = true;
        rng = Math.random;
    } else {
        if(typeof inplace == "function") {
            rng = inplace;
            inplace = true;
        }
    }

    if(inplace) {

        var k = arr.length;

        while(k) {
            var i = Math.floor(rng() * k);
            k--;
            var temp = arr[k];
            arr[k] = arr[i];
            arr[i] = temp;
        }

        return arr;

    } else {

        var result = [];

        for(var i = 0; i < arr.length; ++i) {
            var j = Math.floor(rng() * (i + 1));

            if (j !== i)
                result[i] = result[j];


            result[j] = arr[i];
        }

        return result;
    }
}


//==============================================================================
// Given an integer, returns an array containing its divisors.
//==============================================================================

mishmath.divisors = function(number) {
    let divisors = [];
    const root = Math.floor(Math.sqrt(number));
    for(let i = 1; i <= root; i++) {
        if(number % i === 0) {
            divisors.push(i);
            if(Math.pow(i,2) !== number) {
                divisors.push(number / i);
            }
        }
    }
    return divisors.sort(function(a,b){
        return a - b;
    });
}


//==============================================================================
// Tests an integer for primality, returning true for prime and false for
// composite. Based on is-prime-number by bradleyflood.
//==============================================================================

mishmath.isPrime = function(n) {
    "use strict";

    if(n == 2 || n == 3)
        return true;

    if(n < 2 || n % 1 || n % 2 === 0 || n % 3 === 0)
        return false;

    let root = Math.floor(Math.sqrt(n));

    for(let i = 5; i <= root; i += 6)
        if(n % i === 0 || n % (i + 2) == 0)
            return false;

    return true;
}


//==============================================================================
// Returns the parity or sign of a permutation supplied as an array.
//==============================================================================

mishmath.permutationParity = function(p) {
    "use strict";

    var BRUTE_FORCE_CUTOFF = 32;
    var pool = require("typedarray-pool");

    var n = p.length;

    if(n < BRUTE_FORCE_CUTOFF) { //Use quadratic algorithm for small n

        var sgn = 1;
        for(var i = 0; i < n; ++i) {
            for(var j = 0; j < i; ++j) {
                if(p[i] < p[j]) {
                    sgn = -sgn;
                } else if(p[i] === p[j]) {
                    return 0;
                }
            }
        }
        return sgn;

    } else { //Otherwise use linear time algorithm

        var visited = pool.mallocUint8(n)
        for(var i = 0; i < n; ++i) {
            visited[i] = 0;
        }
        var sgn = 1;
        for(var i = 0; i < n; ++i) {
            if(!visited[i]) {
                var count = 1;
                visited[i] = 1;
                for(var j = p[i]; j !== i; j = p[j]) {
                    if(visited[j]) {
                        pool.freeUint8(visited);
                    return 0;
                    }
                    count += 1;
                    visited[j] = 1;
                }
                if(!(count & 1)) {
                    sgn = -sgn;
                }
            }
        }
        pool.freeUint8(visited);
        return sgn;
    }
}


//==============================================================================
// Given an array of arbitrary elements, arr, returns an array of all
// permutations. If the optional unique argument is true, only unique
// permutations will be returned.
//==============================================================================

mishmath.permutations = function(arr, unique) {
    "use strict";

    if(unique === undefined)
        unique = false;

    let permutations = [];
    let nextPermutation = [];

    function permutate(arr) {
        if (arr.length === 0) {
            permutations.push(nextPermutation.slice());
        }

        for (let i = 0; i < arr.length; i++) {
            arr.push(arr.shift());
            nextPermutation.push(arr[0]);
            permutate(arr.slice(1));
            nextPermutation.pop();
        }
    }

    permutate(arr);

    if(unique)
        return uniq(permutations);

    return permutations;
};


//==============================================================================
// Given an integer, n, returns an array of its prime factors.
//==============================================================================

mishmath.primeFactors = function(n) {
    "use strict";

    var i = 2;
    var res = [];
    var p =0;

    while(n > 1) {
        if(n % i === 0) {
            res.push(i);
            n /= i;
            i = 2;
        } else {
            i++;
        }
    }

    return res;
}


//==============================================================================
// Given polar coordinates r (radius) and t (angle), returns the Cartesian
// coordinates in the form { x: number, y: number }
//==============================================================================

mishmath.polarToCartesian = function(r, t) {
    return { x: r * Math.cos(t), y: r * Math.sin(t) };
}


//==============================================================================
// Given Cartesian coordinates x and y, returns the polar coordinates in the
// form { r: number, t: number } where r is the radial component and t is the
// angular component.
//==============================================================================

mishmath.cartesianToPolar = function(x, y) {
    return { r: Math.sqrt(x*x + y*y), t: Math.atan2(y, x) };
}


//==============================================================================
// Given an angle in degrees, returns its equivalent in radians.
//==============================================================================

mishmath.deg2rad = function(degrees) {
    return degrees * (Math.PI / 180);
}


//==============================================================================
// Given an angle in radians, returns its equivalent in degrees.
//==============================================================================

mishmath.rad2deg = function(radians) {
    return radians * (180 / Math.PI);
}


//==============================================================================
// Given the three-dimensional Cartesian coordinates x, y, z, returns the
// spherical coordinates rho (radial distance), theta (azimuthal angle), and
// phi (polar angle/inclination).
//==============================================================================

mishmath.cartesianToSpherical = function(x, y, z) {
    var rho = Math.sqrt(x*x + y*y + z*z);

    return {
        rho:   rho,
        theta: Math.acos(z / rho),
        phi:   Math.atan(y / x)
    };
}


//==============================================================================
// Given the spherical coordinates rho (radial distance), theta (azimuthal
// angle), and phi (polar angle/inclination), return the three-dimensional Cartesian
// coordinates x, y, z.
//==============================================================================

mishmath.sphericalToCartesian = function(rho, theta, phi) {
    var sinTheta = Math.sin(theta)

    return {
        x: rho * sinTheta * Math.cos(phi),
        y: rho * sinTheta * Math.sin(phi),
        z: rho * Math.cos(theta)
    };
}


//==============================================================================
// Given three-dimensional Cartesian coordinates x, y, z, return the cylindrical
// coordinates rho (radial distance), phi (azimuthal angle in radians),
// and z (axial coordinate).
//==============================================================================

mishmath.cartesianToCylindrical = function(x, y, z) {
    return {
        rho: Math.sqrt(x*x + y*y),
        phi: Math.atan2(y, x),
        z:   z
    };
}


//==============================================================================
// Given cylindrical coordinates rho (radial distance), phi (azimuthal angle in
// radians), and z (axial coordinate), returns the three-dimensional Cartesian
// coordinates x, y, z.
//==============================================================================

mishmath.cylindricalToCartesian = function(rho, phi, z) {
    return {
        x: rho * Math.cos(phi),
        y: rho * Math.sin(phi),
        z: z
    };
}


//==============================================================================
// Converts spherical coordinates rho, theta, and phi into cylindrical coords
// rho, phi, and z.
//==============================================================================

mishmath.sphericalToCylindrical = function(rho, theta, phi) {
    return {
        rho: rho * Math.sin(theta),
        phi: phi,
        z:   rho * Math.cos(theta)
    };
}


//==============================================================================
// Converts cylindrical coordinates rho, phi, and z into spherical coordinates
// rho, theta, and phi.
//==============================================================================

mishmath.cylindricalToSpherical = function(rho, phi, z) {
    return {
        rho:   Math.sqrt(rho*rho + z*z),
        theta: Math.atan(rho/z),
        phi:   phi
    };
}


//==============================================================================
// Given two arrays of coordinates, return the distance between them. Works for
// two or more dimensions. (Euclidean)
//==============================================================================

mishmath.distance = function(a, b) {
    if(a.length == 2)
        return Math.sqrt(Math.pow(a[0] - b[0], 2) + Math.pow(a[1] - b[1], 2));
    else if(a.length == 3)
        return Math.sqrt(Math.pow(a[0] - b[0], 2) + Math.pow(a[1] - b[1], 2) + Math.pow(a[2] - b[2], 2));
    else {
        var sumOfSquares = 0;
        for(var i = 0; i < a.length; i++) {
            sumOfSquares += Math.pow(a[i] - b[i], 2);
        }
        return Math.sqrt(sumOfSquares);
    }
}


//==============================================================================
// Returns a generator which yields successively each combination of M elements
// from a set of N length. The optional mode argument defaults to 'index',
// producing an array of indices. If 'mask' is used, it returns an arrays of
// ones and zeroes.
//==============================================================================

mishmath.combogen = function *(M, N, mode = 'index') {
    var a = new Array(N);
    var c = new Array(M);
    var b = new Array(N);
    var p = new Array(N + 2);
    var x, y, z;

    // init a and b

    for(var i = 0; i < N; i++) {
        a[i] = i;
        if(i < N - M) b[i] = 0;
            else b[i] = 1;
    }

    // init c

    for(i = 0; i < M; i++) {
        c[i] = N - M + i;
    }

    // init p

    for(i = 0; i < p.length; i++) {
        if(i === 0)
            p[i] = N + 1;
        else if(i <= N - M)
            p[i] = 0;
        else if(i <= N)
            p[i] = i - N + M;
        else
            p[i] = -2;
    }

    function twiddle() {
        var i, j, k;
        j = 1;
        while(p[j] <= 0) {
            j++;
        }
        if(p[j - 1] === 0) {
            for(i = j - 1; i !== 1; i--) {
            p[i] = -1;
            }
            p[j] = 0;
            x = z = 0;
            p[1] = 1;
            y = j - 1;
        } else {
            if(j > 1) {
                p[j - 1] = 0;
            }
            do {
                j++;
            } while(p[j] > 0);
            k = j - 1;
            i = j;
            while(p[i] === 0) {
                p[i++] = -1;
            }
            if(p[i] === -1) {
                p[i] = p[k];
                z = p[k] - 1;
                x = i - 1;
                y = k - 1;
                p[k] = -1;
            } else {
                if(i === p[0]) {
                    return 0;
                } else {
                    p[j] = p[i];
                    z = p[i] - 1;
                    p[i] = 0;
                    x = j - 1;
                    y = i - 1;
                }
            }
        }
        return 1;
    }

    if(mode === 'index') {
        yield c.slice();
        while(twiddle()) {
            c[z] = a[x];
            yield c.slice();
        }
    } else if(mode === 'mask') {
        yield b.slice();
        while(twiddle()) {
            b[x] = 1;
            b[y] = 0;
            yield b.slice();
        }
    } else {
        throw new Error('Invalid mode');
    }
};


//==============================================================================
// Given an array of numbers, returns the standard deviation thereof.
//==============================================================================

mishmath.stddev = function(values) {
  "use strict";

  // Get the numbers of data points.
  const length = values.length;

  // Sum the data points.
  const sum = values.reduce((pv, cv) => {
    return pv + cv;
  }, 0);

  // Get the average of the numbers.
  const mean = sum / length;

  // Get the variance (average of the squared differences from the mean)
  const variance = values.reduce((pv, cv) => {
    return pv + Math.pow(cv - mean, 2) / length;
  }, 0);

  // Return the standard Deviation (square root of variance)
  return Math.sqrt(variance);

}


//==============================================================================
// Returns the value of the normal distribution function for a specified value,
// mean, and standard deviation
//==============================================================================

mishmath.normdist = function(x, mean, standardDeviation) {
  "use strict";

  function cdf(x, mean, stdev) {
    return 0.5 * (1 + erf((x - mean) / (Math.sqrt(2 * stdev))));
  }


  function erf(z) {

    const sign = (z >= 0) ? 1 : -1;
    z = Math.abs(z);

    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p  =  0.3275911;

    // A&S formula 7.1.26 with a slice of Horner's
    const t = 1.0/(1.0 + p*z);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-z * z);
    return sign * y;
  }

  return Math.round(100000*cdf(x,mean,standardDeviation))/100000;

}


//==============================================================================
// Given an array of numbers, return their variance.
//==============================================================================

mishmath.variance = function(values) {
    'use strict';

    var mean = mishmath.average(values);

    function sum(a, b) {
        var diff = b - mean;
        return a + (diff * diff);
    }

    return values.reduce(sum, 0) / values.length;
};


//==============================================================================
// Given an array of numbers, return their arithmetic mean.
//==============================================================================

mishmath.average = function(values) {
    'use strict';

    var sum = 0;
    for(var i = 0; i < values.length; i++)
        sum += values[i];

    return sum / values.length;
};


//==============================================================================
// Given two sets of coordinates of any dimensionality, return the Manhattan
// distance between them.
//==============================================================================

mishmath.manhattanDistance = function(a, b) {
    var i = 0,
        ii = a.length,
        d = 0;
    for (; i < ii; i++) {
        d += Math.abs(a[i] - b[i]);
    }
    return d;
};


//==============================================================================
// Given two sets of coordinates of any dimensionality, return the Minkowski
// distance for order p.
//==============================================================================

mishmath.minkowskiDistance = function(a, b, p) {
    var i = 0,
        ii = a.length,
        d = 0;
    for (; i < ii; i++) {
        d += Math.pow(Math.abs(a[i] - b[i]),p);
    }
    return Math.pow(d,(1/p));
};


//==============================================================================
// Given two sets of coordinates of any dimensionality, return the Chebyshev
// distance between them.
//==============================================================================

mishmath.chebyshevDistance = function(a, b) {
    var ii = a.length,
        max = 0,
        aux = 0;
    for (var i = 0; i < ii ; i++) {
        aux = Math.abs(a[i] - b[i]);
        if (max < aux) {
            max = aux;
        }
    }
    return max;
};


//==============================================================================
// Calculates the distance between two sets of decimal latitude and longitude
// coordinates on the surface of Earth using the Haversine formula. By default,
// the result is in meters. To use different units (or to calculate distances on
// another planet entirely), supply an explicit value for planetRadius in the
// units desired.
//==============================================================================

mishmath.haversineDistance = function(lat1, lon1, lat2, lon2, planetRadius = 6371000) {

  var latitudeDifference = this.deg2rad(lat2 - lat1);
  var logitudeDifference = this.deg2rad(lon2 - lon1);

  var n =
    Math.sin(latitudeDifference / 2) * Math.sin(latitudeDifference / 2) +
    Math.cos(this.deg2rad(lat1)) * Math.cos(this.deg2rad(lat2)) *
    Math.sin(logitudeDifference / 2) * Math.sin(logitudeDifference / 2);

  var distance = 2 * Math.atan2(Math.sqrt(n), Math.sqrt(1 - n));

  return planetRadius * distance;
};


//==============================================================================
//
// Ported from Mukesh Prasad's public domain code:
//   http://tog.acm.org/resources/GraphicsGems/gemsii/xlines.c
//
//  This function computes whether two line segments,
//  respectively joining the input points (x1,y1) -- (x2,y2)
//  and the input points (x3,y3) -- (x4,y4) intersect.
//  If the lines intersect, the return value is an array
//  containing coordinates of the point of intersection.
//
//  Params
//       x1, y1,  x2, y2   Coordinates of endpoints of one segment.
//       x3, y3,  x4, y4   Coordinates of endpoints of other segment.
//
//  Also Accepts:
//   4 objects with the minimal object structure { x: .., y: ..}
//   4 arrays where [0] is x and [1] is y
//
//  The value returned by the function is one of:
//
//       undefined - no intersection
//       array     - intersection
//       true      - colinear
//==============================================================================

mishmath.segmentsIntersect = function (x1, y1, x2, y2, x3, y3, x4, y4) {

    if (arguments.length === 4) {
        var p1 = x1;
        var p2 = y1;
        var p3 = x2;
        var p4 = y2;

        // assume array [x, y]

        if (p1.length && p1.length === 2) {
            x1 = p1[0];
            y1 = p1[1];
            x2 = p2[0];
            y2 = p2[1];
            x3 = p3[0];
            y3 = p3[1];
            x4 = p4[0];
            y4 = p4[1];

        // assume object with obj.x and obj.y

        } else {
            x1 = p1.x;
            y1 = p1.y;
            x2 = p2.x;
            y2 = p2.y;
            x3 = p3.x;
            y3 = p3.y;
            x4 = p4.x;
            y4 = p4.y;
        }
    }


    var a1, a2, b1, b2, c1, c2; // Coefficients of line eqns.
    var r1, r2, r3, r4;         // 'Sign' values
    var denom, offset;          // Intermediate values
    var x, y;                   // Intermediate return values

    // Compute a1, b1, c1, where line joining points 1 and 2
    // is "a1 x  +  b1 y  +  c1  =  0".

    a1 = y2 - y1;
    b1 = x1 - x2;
    c1 = x2 * y1 - x1 * y2;

    // Compute r3 and r4.

    r3 = a1 * x3 + b1 * y3 + c1;
    r4 = a1 * x4 + b1 * y4 + c1;

    // Check signs of r3 and r4.  If both point 3 and point 4 lie on
    // same side of line 1, the line segments do not intersect.

    if ( r3 !== 0 && r4 !== 0 && ((r3 >= 0 && r4 >= 0) || (r3 < 0 && r4 < 0))) {
        return; // no intersection
    }

    // Compute a2, b2, c2

    a2 = y4 - y3;
    b2 = x3 - x4;
    c2 = x4 * y3 - x3 * y4;

    // Compute r1 and r2

    r1 = a2 * x1 + b2 * y1 + c2;
    r2 = a2 * x2 + b2 * y2 + c2;

    // Check signs of r1 and r2.  If both point 1 and point 2 lie
    // on same side of second line segment, the line segments do
    // not intersect.

    if (r1 !== 0 && r2 !== 0 && ((r1 >= 0 && r2 >= 0) || (r1 < 0 && r2 < 0))) {
        return; // no intersections
    }

    // Line segments intersect: compute intersection point.

    denom = a1 * b2 - a2 * b1;

    if ( denom === 0 ) {
        return true;
    }

    offset = denom < 0 ? - denom / 2 : denom / 2;

    x = b1 * c2 - b2 * c1;
    y = a2 * c1 - a1 * c2;

    return [
        ( x < 0 ? x : x ) / denom,
        ( y < 0 ? y : y ) / denom,
    ];
};


//==============================================================================
// Returns the shortest distance from point (x, y) to line segment (x1, y1) -
// (x2, y2). By https://stackoverflow.com/users/368954/joshua.
//==============================================================================

mishmath.pointToLine = function(x, y, x1, y1, x2, y2) {

    var A = x - x1;
    var B = y - y1;
    var C = x2 - x1;
    var D = y2 - y1;

    var dot = A * C + B * D;
    var len_sq = C * C + D * D;
    var param = -1;
    if (len_sq != 0) //in case of 0 length line
        param = dot / len_sq;

    var xx, yy;

    if(param < 0) {
        xx = x1;
        yy = y1;
    } else if (param > 1) {
        xx = x2;
        yy = y2;
    } else {
        xx = x1 + param * C;
        yy = y1 + param * D;
    }

    var dx = x - xx;
    var dy = y - yy;
    return Math.sqrt(dx * dx + dy * dy);
}


//==============================================================================
// Returns a boolean indicating whether point is inside vs, where point is a
// two-element array and vs is an array of two-element arrays or a flat array.
//
// ray-casting algorithm based on
// https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html

mishmath.pointInPolygon = function(point, vs, start, end) {
    if (vs.length > 0 && Array.isArray(vs[0])) {
        return pointInPolygonNested(point, vs, start, end);
    } else {
        return pointInPolygonFlat(point, vs, start, end);
    }
}

//------------------------------------------------------------------------------

function pointInPolygonFlat (point, vs, start, end) {
    var x = point[0], y = point[1];
    var inside = false;
    if (start === undefined) start = 0;
    if (end === undefined) end = vs.length;
    var len = (end-start)/2;
    for (var i = 0, j = len - 1; i < len; j = i++) {
        var xi = vs[start+i*2+0], yi = vs[start+i*2+1];
        var xj = vs[start+j*2+0], yj = vs[start+j*2+1];
        var intersect = ((yi > y) !== (yj > y))
            && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    return inside;
};

//------------------------------------------------------------------------------

function pointInPolygonNested (point, vs, start, end) {
    var x = point[0], y = point[1];
    var inside = false;
    if (start === undefined) start = 0;
    if (end === undefined) end = vs.length;
    var len = end - start;
    for (var i = 0, j = len - 1; i < len; j = i++) {
        var xi = vs[i+start][0], yi = vs[i+start][1];
        var xj = vs[j+start][0], yj = vs[j+start][1];
        var intersect = ((yi > y) !== (yj > y))
            && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    return inside;
};


//==============================================================================
// Given an array of n-element point arrays, return the centroid.
//==============================================================================

mishmath.centroid = function(points) {
  if (typeof points === 'undefined') {
    return [];
  }

  let dimensions = points[0].length;
  let accumulation = points.reduce((acc, point) => {
    point.forEach((dimension, idx) => {
      acc[idx] += dimension;
    });

    return acc;
  }, Array(dimensions).fill(0));

  return accumulation.map(dimension => dimension / points.length);
}


//==============================================================================
// Returns the area of a polygon specified as an array of vertices.
//==============================================================================

mishmath.polygonArea = function(polygon) {
    return Math.abs(_polygonArea(polygon));
}


//==============================================================================
// Given a polygon specified as an array of vertices, returns boolean true if
// they are in clockwise order or false if they are in counterclockwise order.
//==============================================================================

mishmath.polygonIsClockwise = function(polygon) {
    return _polygonArea(polygon) < 0 ? true : false;
}

//------------------------------------------------------------------------------

function _polygonArea(a) {
  var e0 = [0, 0];
  var e1 = [0, 0];
  var area = 0;
  var first = a[0];

  var l = a.length;
  for (var i=2; i<l; i++) {
    var p = a[i-1];
    var c = a[i];
    e0[0] = first[0] - c[0];
    e0[1] = first[1] - c[1];
    e1[0] = first[0] - p[0];
    e1[1] = first[1] - p[1];

    area += (e0[0] * e1[1]) - (e0[1] * e1[0]);
  }
  return area/2;
}


//==============================================================================
// Given two polygons in the form of two arrays of x,y coordinates, returns a
// boolean indicating whether they overlap.

mishmath.polygonOverlap = function(polyA, polyB) {

    // First, check for points in polygons -------------------------------------

    for(var p = 0; p < polyA.length; p++)
        if(mishmath.pointInPolygon(polyA[p], polyB))
            return true;

    for(var p = 0; p < polyB.length; p++)
        if(mishmath.pointInPolygon(polyB[p], polyA))
            return true;

    // Second, check for line intersections ------------------------------------

    for(var p1 = 0; p1 < polyA.length; p1++) {
        for(var p2 = 0; p2 < polyB.length; p2++) {
            var res = mishmath.segmentsIntersect(
                polyA[p1][0], polyA[p1][1],
                polyA[(p1+1) % polyA.length][0], polyA[(p1+1) % polyA.length][1],
                polyB[p2][0], polyB[p2][1],
                polyB[(p2+1) % polyB.length][0], polyB[(p2+1) % polyB.length][1]
            );
            if(Array.isArray(res))
                return true;
        }
    }

    return false;

}


//==============================================================================

module.exports = mishmath;

