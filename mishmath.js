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








module.exports = mishmath;

