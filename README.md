# mishmath v0.0.1

**A collection of miscellaneous math routines for Node.js, mostly culled from other FOSS modules**

<a name="introduction"></a>
## Introduction

There are a lot of good but abandoned npm modules for various useful math and 
logic functions out there, and I find myself using them -- or at least thinking 
I'll have a use for them -- frequently enough that I've started lumping them 
together into a generic module for my own use. Mishmath is an attempt to clean 
up, organize, and maintain that module for others to use.

You might wonder if one monolithic module is preferable to dozens of independent 
modules that you can use as needed. It is preferable in terms of convenience, 
but not in terms of space efficiency. For that reason, I have endeavored to 
structure the code in such a way that it is easy to pull out what you need and 
stuff it into a custom module if that is what you need. Additionally, I have
taken pains to document the original source modules (and credit their authors)
so it is easy to go back to the sources.

<a name="what"></a>
## What's in Mishmath?

As the terrible pun suggests, Mishmath is a heterogenous collection of routines 
I found useful. It is both general and highly idiosyncratic, aiming neither to 
specialize in any one area or to be comprehensive. A lot of what I do 
professionally is related to combinatorics and geometry, so there's a lot of 
that, plus a bunch of obscure odds and ends. It's also continuously growing 
based on my current needs, so there's no telling what will be next.

Suggestions and contributions are always welcome.


<a name="functions"></a>
## Functions

**isPrime(n)**

Tests `n` for primality, returning `true` if prime or `false` if composite.



<a name="Credits"></a>
## Credits

I cannot be thankful enough to the many people who have built useful, quality
code and given it to the larger community of Node developers under permissive
open source licenses. 

**isPrime** is derived from [Bradley Flood](https://www.npmjs.com/~bradleyflood)'s
[is-prime-number](https://www.npmjs.com/package/is-prime-number).


<a name="license"></a>
## License

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

<a name="todo"></a>
## Todo

* Everything
* Solicit suggestions from users.

<a name="changelog"></a>
## Changelog

0.0.1: Initial release.
