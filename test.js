#!/usr/bin/env node

const mm = require("./mishmath.js");

// divisors ====================================================================

function testDivisors() {
    for(var n = 0; n < 26; n++) {
        var div = mm.divisors(n);
        console.log("divisors -- divisors of " + n + ": ", div);
    }
}


// isPrime =====================================================================

function testIsPrime() {
    var primes = [ ];
    var n = 1;
    while(primes.length < 100) {
        if(mm.isPrime(n))
            primes.push(n);
        n++;
    }
    console.log("isPrime -- first hundred primes: " + primes.join(", "));
}


// permutations ================================================================

function testPermutations() {
    console.log("testPermutations -- Permutations of ['a', 'b', 'c']:\n", mm.permutations(['a', 'b', 'c']));
    console.log("testPermutations -- Permutations of [1, 2, 3]:\n", mm.permutations([1, 2, 3]));
    console.log("testPermutations -- Unique permutations of ['a', 'a', 'b']:\n", mm.permutations(['a', 'a', 'b'], true));
}


// permutationParity ===========================================================

function testPermutationParity() {
    var arr = mm.permutations(['a', 'b', 'c']);
    for(var i = 0; i < arr.length; i++) {
        console.log("permutationParity -- Parity of " + arr[i] + " = " + mm.permutationParity(arr[i]));
    }
}


// primeFactors ================================================================

function testPrimeFactors() {
    for(var n = 0; n < 26; n++) {
        var pf = mm.primeFactors(n);
        console.log("primeFactors -- prime factors of " + n + ": ", pf);
    }
}


// cartesianToPolar/polarToCartesian============================================

function testCartesianToPolar() {
    for(var x = 10; x > -11; x -= 2) {
        for(var y = -10; y < 11; y += 1.5) {
            var p = mm.cartesianToPolar(x, y);
            console.log("cartesianToPolar: x,y " + x + "," + y + " = r,t " + p.r + "," + p.t);
            var c = mm.polarToCartesian(p.r, p.t);
            console.log("polar2Cartesian: r,t " + p.r + "," + p.t + " = x,y " + c.x + "," + c.y );
        }
    }
}

// deg2rad/rad2deg =============================================================

function testDeg2rad() {
    for(var degrees = 0; degrees < 361; degrees += 10) {
        var radians = mm.deg2rad(degrees);
        console.log("deg2rad: " + degrees + " degrees = " + radians + " radians.");
        var redeg = mm.rad2deg(radians);
        console.log("rad2deg: " + radians + " radians = " + redeg + " degrees.");
    }
}

// cartesianToSpherical/sphericalToCartesian ===================================

function testCartesianSpherical() {
    var cc = [
        [10, 1, 8 ],
        [10, 2, 7 ],
        [10, 3, 6 ],
        [10, 4, 5 ],
        [-10, 5, 4 ],
        [-10, 6, 3 ],
        [-10, 7, 2 ],
        [-10, 8, 1 ]
    ];

    for(var i = 0; i < cc.length; i++) {
        var sph = mm.cartesianToSpherical(cc[i][0], cc[i][1], cc[i][2]);
        console.log("cartesianToSpherical: " + cc[i].join(",") + " --> " + [sph.rho, sph.theta, sph.phi] + ".");
        var cart = mm.sphericalToCartesian(sph.rho, sph.theta, sph.phi);
        console.log("sphericalToCartesian: " + [sph.rho, sph.theta, sph.phi] + " --> " + [cart.x, cart.y, cart.z] + ".");
    }
}

// cartesianToCylindrical/cylindricalToCartesian ===============================

function testCartesianCylindrical() {
    var cc = [
        [10, 1, 8 ],
        [10, 2, 7 ],
        [10, 3, 6 ],
        [10, 4, 5 ],
        [-10, 5, 4 ],
        [-10, 6, 3 ],
        [-10, 7, 2 ],
        [-10, 8, 1 ]
    ];

    for(var i = 0; i < cc.length; i++) {
        var cyl = mm.cartesianToCylindrical(cc[i][0], cc[i][1], cc[i][2]);
        console.log("cartesianToCylindrical: " + cc[i].join(",") + " --> " + [cyl.rho, cyl.phi, cyl.z] + ".");
        var cart = mm.cylindricalToCartesian(cyl.rho, cyl.phi, cyl.z);
        console.log("cylindricalToCartesian: " + [cyl.rho, cyl.phi, cyl.z] + " --> " + [cart.x, cart.y, cart.z] + ".");
    }
}

// sphericalToCylindrical/cylindricalToSpherical ===============================

function testSphericalCylindrical() {
    var orig = { rho: 12.845, theta: 0.898, phi: 0.0996 };
    var cyl = mm.sphericalToCylindrical(orig.rho, orig.theta, orig.phi);
    var sph = mm.cylindricalToSpherical(cyl.rho, cyl.phi, cyl.z);
    console.log("spherical to cylindrical and back: ", orig, cyl, sph);
}


// fisherYatesShuffle =========================================================


function testFisherYatesShuffle() {
    var arr1 = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ];
    var arr2;

    for(var n = 0; n < 10; n++) {
        mm.fisherYatesShuffle(arr1);
        console.log("fisherYatesShuffle/inplace: ", arr1);
    }

    for(var n = 0; n < 10; n++) {
        arr2 = mm.fisherYatesShuffle(arr1);
        console.log("fisherYatesShuffle/copy: ", arr2);
    }
}

// distance ====================================================================

function testDistance() {
    console.log("Diagonal of a unit square:    " + mm.distance([0, 0], [1, 1]));
    console.log("Diagonal of a unit cube:      " + mm.distance([0,0,0], [1,1,1]));
    console.log("Diagonal of a unit hypercube: " + mm.distance([0,0,0,0], [1,1,1,1]));
}

// combogen ====================================================================

function testCombogen() {
    var gen = mm.combogen(2, 4);
    console.log([...gen]);
    gen = mm.combogen(2, 4, 'mask');
    console.log([...gen]);
}

// stddev ======================================================================

function testStddev() {
    console.log(mm.stddev([600,470,170,430,300])); // 147.32277488562318
    console.log(mm.stddev([13,23,12,44,55])); // 17.21162397916013
}

// normdist =====================================================================

function testNormdist() {
    console.log(mm.normdist(2,0,1)); // 0.97725
    console.log(mm.normdist(0,0,1)); // 0.5
}

// average ======================================================================

function testAverage() {
    console.log(mm.average([100, 100, 50, 50]));
}

// variance =====================================================================

function testVariance() {
    console.log(mm.variance([600, 470, 170, 430, 300]));
}

// manhattanDistance ============================================================

function testManhattanDistance() {
    for(i = 0; i < 6; i++) {
        console.log("Manhattan distance from origin to " + i + ", " + i + ": " + mm.manhattanDistance([0,0], [i,i]));
    }
}

// minkowskiDistance ============================================================

function testMinkowskiDistance() {
    for(i = 0; i < 6; i++) {
        console.log("Minkowski distance of order 1 from origin to " + i + ", " + i + ": " + mm.minkowskiDistance([0,0], [i,i], 1));
        console.log("Minkowski distance of order 2 from origin to " + i + ", " + i + ": " + mm.minkowskiDistance([0,0], [i,i], 2));
    }

}

// chebyshevDistance ============================================================

function testChebyshevDistance() {
    for(i = 0; i < 6; i++) {
        console.log("Chebyshev distance from origin to " + i + ", " + (i - 1) + ": " + mm.chebyshevDistance([0,0], [i,i-1]));
    }
}


testDivisors();
testIsPrime();
testPermutations();
testPermutationParity();
testPrimeFactors();
testCartesianToPolar();
testDeg2rad();
testCartesianSpherical();
testCartesianCylindrical();
testSphericalCylindrical();
testFisherYatesShuffle();
testDistance();
testCombogen();
testStddev();
testNormdist();
testAverage();
testVariance();
testManhattanDistance();
testMinkowskiDistance();
testChebyshevDistance();