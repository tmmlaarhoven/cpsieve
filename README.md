# The CPSieve algorithm

Included in this repository is a proof-of-concept implementation of the CPSieve algorithm for solving the shortest vector problem (SVP) on lattices, introduced in the paper:

 - Anja Becker, Thijs Laarhoven, **Efficient (ideal) lattice sieving using cross-polytope LSH**, *Cryptology ePrint Archive, Report 2015/823*.

This code is released under the MIT license: you can basically do anything you want with/to the code, as long as you give proper attribution back to us, and as long as you do not hold us responsible for anything that happens to your machine as a result of running this code.

### Arbitrary lattices

For solving the shortest vector problem on arbitrary lattices, use the file `cpsieve.cpp`. This file is currently set to verbose mode, outputting new records each time it finds a new shorter vector in the input lattice. The input filename is assumed to follow the format "dim40sd0.txt", as is also the case for the [SVP challenge], but this format can be adjusted at the beginning of the main algorithm loop. The algorithm further keeps track of various statistics, which are stored in the file `cpsieve-results.txt`. This implementation should solve SVP correctly for *any* basis of a lattice, if it is given enough time and space to perform its operations. (Note that things may start to break down if the coordinates or norms of vectors become too large, which is not an issue for the SVP challenge lattices.)

### Ideal lattices

For solving SVP efficiently on ideal lattices, there are two faster versions of the above implementation, as described in the paper; one for ideals with modulus `X^n + 1`, and one for ideals with modulus `X^n - 1`.

##### Power-of-2 cyclotomic lattices (`X^n + 1`)

For lattices where the underlying ideal is taken modulo `X^n + 1` with `n = 2^k` a power of 2, as commonly used in lattice-based cryptography, the implementation given in `icpsieve.cpp` achieves better results. This implementation exploits the additional structure to achieve a linear speed-up and a linear decrease in the space complexity. The input basis file is assumed to be of the form `idealindex64dim32sd0.txt` where `32 = (64 / 2)` is the dimension of the lattice and `0` is the random seed. These challenge bases can for instance be downloaded from the [ideal SVP challenge].

##### NTRU lattices (`X^n - 1`)

For NTRU lattices, similar time and memory improvements can be achieved using techniques described in the paper. For NTRU lattices, where the underlying polynomial ring is given modulo a prime `q`, and modulo a polynomial `X^m - 1` (where `n = 2 m` is the dimension of the block-cyclic NTRU lattice), one can use the implementation given in `ìcpsieve-ntru.cpp`. As described in the paper, this also leads to considerable improvements, although the improvements are a bit smaller than for the power-of-2 cyclotomic lattices. The input basis file is assumed to be of the form `NTRU31seed0.txt` where the prime `m = 31` is half the dimension of the lattice and `0` is the random seed. NTRU lattices can be generated using e.g. [Sage].

For any further questions or comments, feel free to contact us at mail@thijs.com or anja.becker@epfl.ch.

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)
  
  [SVP challenge]: http://latticechallenge.org/svp-challenge/index.php
  [ideal SVP challenge]: http://latticechallenge.org/svp-challenge/index.php
  [Sage]: http://www.sagemath.org/

