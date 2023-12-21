# mishmath v0.0.11

![mishmath title](img/mishmath.png)

**A collection of miscellaneous math routines for Node.js, mostly culled from other FOSS modules**

**NEW in 0.0.11**: Polygon overlap test.

## Table of Contents

* [Introduction](#introduction)
* [What's in Mishmath?](#what)
* [Functions](#functions)
* [Credits](#credits)
* [License](#license)
* [Todo](#todo)
* [Changelog](#changelog)


<a name="introduction"></a>
## Introduction

There are a lot of good but abandoned npm modules for various useful math and
logic functions out there, and I find myself using them -- or at least thinking
I'll have a use for them -- frequently enough that I've started lumping them
together into a generic module for my own use. Mishmath is an attempt to clean
up, organize, and maintain a public module for others to use.

You might wonder if one monolithic module is preferable to dozens of independent
modules that you can use as needed. It is preferable in terms of convenience,
but not in terms of space efficiency. For that reason, I have endeavored to
structure the code in such a way that it is easy to pull out what you need and
stuff it into a custom module. Additionally, I have taken pains to document the
original source modules (and credit their awesome authors) so it is easy to go
back to the sources.

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

**average(ary)**

Given an array of numbers, returns their average or arithmetic mean.

---

**cartesianToPolar(x, y)**

Given Cartesian coordinates `x` and `y`, returns an object of the form
`{r: num, t: num} where `r` is the radial/rho value and `t` is
the angular/theta value (in radians).

---

**cartesianToSpherical(x, y, z)**

Given 3D Cartesian coordinates `x`, `y`, `z`, returns an object of the form
`{rho: num, theta: num, phi: num}` where `rho` is the radial distance, `theta`
is the azimuthal angle in radians, and `phi` is polar angle (inclination) in
radians.

---

**centroid(points)**

Given an array of `points`, each point being an array of arbitrary (but equal)
length, return their centroid.

---

**chebyshevDistance(a, b)**

Given two coordinate vectors of any dimensionality, return the Chebyshev
distance between them.

---

**combogen(m, n, mode)**

This little gem returns a generator function which yields successively each
combination of `m` elements from a set of `n` members. If the optional `mode`
argument is `'index'` (the default), it returns `m`-length arrays of indices. If
`mode` is `'mask'`, it returns `n`-length arrays of ones and zeroes.

```javascript
var gen = combogen(2, 4);
console.log([...gen]); // [ [ 3, 2 ], [ 0, 2 ], [ 1, 2 ], [ 1, 2 ], [ 0, 2 ], [ 0, 1 ] ]

gen = combogen(2, 4, 'mask');
console.log([...gen]); // [ [ 0, 0, 1, 1 ],[ 1, 0, 0, 1 ],[ 0, 1, 0, 1 ],[ 0, 1, 1, 0 ],[ 1, 0, 1, 0 ],[ 1, 1, 0, 0 ] ]
```

---

**cylindricalToSpherical(rho, phi, z)**

Given cylindrical coordinates `rho` (radial distance), `phi` (azimuth), and `z`
(height), returns an object of the form `{rho: num, theta: num, phi: num}` where
`rho` is the radial distance, `theta` is the azimuthal angle in radians, and
`phi` is polar angle (inclination) in radians.

---

**deg2rad(degrees)**

Given an angle in `degrees`, returns its equivalent in radians.

---

**distance(a, b)**

Returns the Euclidean distance between two points, `a` and `b`, which are both
represented as arrays of coordinates. Works for any arbitrary dimension of 2 or
higher.

---

**divisors(n)**

Given an integer, `n`, return an array containing all of its divisors.

---

**fisherYatesShuffle(arr, inplace, rng)**

Performs the fast Fisher-Yates randomization of an array, `arr`. If the optional
`inplace` argument is `true` (the default), `arr` is sorted in place; otherwise,
a new array is created. The optional `rng` argument can be used to supply a
random number generation callback to replace the default `Math.random`. Returns
the shuffled array.

---

**haversineDistance(lat1, lon1, lat2, lon2, planetRadius = 6371000)**

Calculates the distance between two sets of decimal latitude and longitude
coordinates on the surface of Earth using the Haversine formula. By default, the
result is in meters. To use different units (or to calculate distances on
another planet entirely), supply an explicit value for `planetRadius` in the
units desired.

---

**isPrime(n)**

Tests `n` for primality, returning `true` if prime or `false` if composite.

---

**manhattanDistance(a, b)**

Given two coordinate vectors of any dimensionality, return the Manhattan
distance between them.

---

**minkowskiDistance(a, b, p)**

Given two coordinate vectors of any dimensionality, return the Minkowski
distance between them for order `p`. For `p == 1`, this is equivalent to the
Manhattan distance, and `p == 2` is equivalent to the Euclidean distance.

---

**normdist(value, mean, stddev)**

Given a value, mean, and standard deviation, returns the value of the normal
distribution.

---

**permutationParity(arr)**

Given a permutation in the form of an array, returns its [parity or
sign](https://en.wikipedia.org/wiki/Parity_of_a_permutation). The parity is
represented by `1` if the permutation is odd, `-1` if it is even, and `0` if
`arr` is not a permutation

---

**permutations(arr, unique = false)**

Given an array of arbitrary elements, `arr`, returns an array of all
permutations. If the optional `unique` argument is `true`, only unique
permutations will be returned.

---

**pointInPolygon(point, polygon, start = 0, end = polygon.length)**

Determines whether `point` is inside `polygon` using a ray-casting algorithm.
The `point` is given as a two-element array, while `polygon` may be either
an array of two-element points or a flat array of coordinates.

The optional `start` and `end` coordinates allow you to specify a range of
coordinates in `polygon` to use for the actual polygon in cases where many
polygons are specified in the actual array.

---

**pointToLine(x, y, x1, y1, x2, y2)**

Returns the shortest distance from point (`x`, `y`) to line segment (`x1`, `y1`) -
(`x2`, `y2`).

---

**polarToCartesian(r, t)**

Given polar coordinates `r` (radial/rho) and `t` (angular/theta), the latter
in radians, returns an object of the form `{x: num, y: num}`.

---

**polygonArea(p)**

Returns the area of a polygon specified as an array of `[x,y]` vertices.

---

**polygonIsClockwise(p)**

Given a polygon specified as an array of `[x,y]` vertices, returns boolean
`true` if the vertices are in clockwise order, `false` otherwise.

---

**polygonOverlap(poly1, poly2)**

Given two polygons specified as arrays of `[x,y]` vertices, returns boolean
`true` if the polygons overlap, `false` otherwise. N.b.: Colinear polygons
are not considered overlaps here.

---

**primeFactors(n)**

Returns an array containing the prime factors of `n`.

---

**rad2deg(radians)**

Given an angle in `radians`, returns its equivalent in degrees.

---

**segmentsIntersect(x1, y1, x2, y2, x3, y3, x4, y4)**

Given two line segments (`x1`, `y1`) -- (`x2`, `y2`) and (`x3`, `y3`) -- (`x4`, `y4`),
returns an array with the coordinates of their intersection if they in
fact intersect. If they are colinear, boolean `true` is returned, and
if they neither intersect nor are colinear, `undefined` is returned.

The arguments may also be passed as four objects of the form `{ x: ..., y: ... }`
or four `[x, y]` arrays.

---

**sphericalToCartesian(rho, theta, phi)**

Given spherical coordinates `rho` (radial distance), `theta` (azimuthal angle in
radians), and `phi` (polar angle or inclination in radians), returns an object
of the form `{x: num, y: num, z: num}`.

---

**sphericalToCylindrical(rho, theta, phi)**

Given Given spherical coordinates `rho` (radial distance), `theta` (azimuthal angle in
radians), and `phi` (polar angle or inclination in radians), returns an object
of the form `{rho: num, phi: num, z: num}` where `rho` is the radial distance, `phi`
is the azimuth, and `z` is the height.

---

**stddev(ary)**

Takes an array of numbers and returns their standard deviation.

---

**variance(ary)**

Given an array of numbers, returns the array's variance.

<a name="Credits"></a>
## Credits

I cannot be thankful enough to the many people who have built useful, quality
code and given it to the larger community of Node developers under permissive
open source licenses.

Even if you don't care about the original version of the functions in `mishmath`,
the developers' profile pages are worth looking at. Many of them are prolific
module authors, and all of them are quite good.

**combogen** is based on the very nicely coded [`ml-combinations`](https://www.npmjs.com/package/ml-combinations)
package.

**centroid** is based on [Joe D'Alessandro](https://www.npmjs.com/~xupit3r)'s
elegant and compact [centroider](https://www.npmjs.com/package/centroider).

**divisors** is based on [janjarfalk](https://www.npmjs.com/~janjarfalk)'s
[`get-divisors`](https://www.npmjs.com/package/get-divisors).

**fisherYatesShuffle** is derived from [dcousens](https://www.npmjs.com/~dcousens)'
[`fisher-yates`](https://www.npmjs.com/package/fisher-yates).

**haversineDistance** was lifted from the tidy code in
[seabass](https://www.npmjs.com/~seabass)'s
[s-haversine](https://www.npmjs.com/package/s-haversine) package.

**isPrime** is derived from [Bradley Flood](https://www.npmjs.com/~bradleyflood)'s
[`is-prime-number`](https://www.npmjs.com/package/is-prime-number). The original
had a subtle bug in it that miscategorized 2 and 3 as composite numbers, which
has now been fixed.

**manhattanDistance**, **minkowskiDistance**, and **chebyshevDistance** are from
[Michaël Zasso](https://github.com/targos) and [Miguel Angel Asencio Hurtado](https://github.com/maasencioh)'s
[`ml-distance`](https://www.npmjs.com/package/ml-distance), which
is a smorgasbord of distance and similarity metrics.

**normdist** is from [janjarfalk](https://www.npmjs.com/~janjarfalk)'s
[`get-normal-distribution`](https://www.npmjs.com/package/get-normal-distribution).

**permutationParity** was pulled out of [Mikola Lysenko](https://www.npmjs.com/~mikolalysenko)'s
[`permutation-parity`](https://www.npmjs.com/package/permutation-parity). Mikola has
a mountain of other awesome math/graphics modules that are well worth checking out.

**permutations** is based on [janjarfalk](https://www.npmjs.com/~janjarfalk)'s
[`get-permutations`](https://www.npmjs.com/package/get-permutations). The `uniq`
function used for the unique permutations option came out of janjarfalk's
[`get-unique-permutations`](https://www.npmjs.com/package/get-unique-permutations),
where it is wrapped around the original `get-permutations` module.

**pointInPolygon** is from [nopersonsmodules](https://www.npmjs.com/~nopersonsmodules)'s
elegant [point-in-polygon](https://www.npmjs.com/package/point-in-polygon), which is
in turn based this [ray-casting algorithm](https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html).

**pointToLine** is one of the few I came across on StackOverflow, which does contain
some gems if you're willing to sort through a few megatons of overconfident people
incorrecting each other. This gem is by [Joshua Perina](http://joshua.perina.com/).

**polygonArea** and **polygonIsClockwise** are derived from
[tmpvar](https://www.npmjs.com/~tmpvar)'s [2d-polygon-area](https://www.npmjs.com/package/2d-polygon-area).

**primeFactors** is derived from [janjarfalk](https://www.npmjs.com/~janjarfalk)'s
[`get-prime-factors`](https://www.npmjs.com/package/get-prime-factors).

**segmentsIntersect** is based on [tmpvar](https://www.npmjs.com/~tmpvar)'s
[segseg](https://www.npmjs.com/package/segseg), which in turn is based on
Mukesh Prasad's public domain code in _Graphics Gems_, plus some ease-of-use
additions (by tmpvar, not me).

**stddev** came from [janjarfalk](https://www.npmjs.com/~janjarfalk)'s
[`get-standard-deviation`](https://www.npmjs.com/package/get-standard-deviation).

**variance** was taken from [bytespider](https://www.npmjs.com/~bytespider)'s
[`variance`](https://www.npmjs.com/package/variance) package.

<a name="license"></a>
## License

Copyright 2019-2023 Eric O'Dell, original authors, and subsequent contributors

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

* Polygon self-intersection
* Spline interpolation
* Continue to cherry-pick ml-distance?
* Generatorics
* Wrapper around combogen to produce explicit array
* Median
* Linear regression
* Lehmer codes

<!--

npm packages to plunder:

disi
foreach-combination
furlong
generatorics
get-average-if
get-matrix-sums
historical-permutations
js-combinations
ml-distance
ordered-char-combinations
permutron
remap-value
typedarray-pool

npm packages that are too popular and/or not suited for mishmath, but are still awesome.

curve-interpolator
wuzzy

-->

<a name="changelog"></a>
## Changelog

0.0.9:
* `polygonArea`
* `polygonIsClockwise`

0.0.8:

* Removed temporary dependencies.
* `pointInPolygon`
* `centroid`

0.0.7:

* `pointToLine`

0.0.6:

* `haversineDistance`
* `segmentsIntersect`

0.0.5:

* `manhattanDistance`
* `minkowskiDistance`

0.0.4: Weekly update (2019-06-06), including:

* `average`
* `combogen`
* `distance`
* `normdist`
* `variance`

0.0.3: Weekly update (2019-05-28), including

* `cartesianToCylindrical`
* `cartesianToPolar`
* `cartesianToSpherical`
* `cylindricalToCartesian`
* `cylindricalToSpherical`
* `deg2rad`
* `polarToCartesian`
* `rad2deg`
* `sphericalToCartesian`
* `sphericalToCylindrical`

0.0.2: Fat fingers, ignore.

0.0.1: Initial release, including

* `divisors`
* `fisherYatesShuffle`
* `isPrime`
* `permutationParity`
* `permutations` plus unique option
* `primeFactors`

