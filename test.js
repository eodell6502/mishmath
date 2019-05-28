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




