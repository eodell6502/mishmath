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

//testDivisors();
testFisherYatesShuffle();
//testIsPrime();
//testPermutations();
//testPermutationParity();
//testPrimeFactors();

