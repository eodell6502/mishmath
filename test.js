const mm = require("./mishmath.js");

// isPrime =====================================================================

function testIsPrime() {
    var primes = [ ];
    var n = 2;
    while(primes.length < 100) {
        if(mm.isPrime(n))
            primes.push(n);
        n++;
    }
    console.log("isPrime - first hundred primes: " + primes.join(", "));
}



testIsPrime();
