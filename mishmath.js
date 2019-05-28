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


module.exports = mishmath;

