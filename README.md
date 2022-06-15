# fastpoly

Proof of concept implementation of fast polynomial multiplication using Fourier transform

Originally written 12/13/2020

The program implements fast polynomial multiplication as well as the standard "slow" O(n^2) algorithm.
It generates random polynomials of a given degree and compares the time taken.
For small degree polynomials, the "slow" multiplication method is still quicker.

Inspired by (this video by Reducible)[https://www.youtube.com/watch?v=h7apO7q16V0].

## Sample Output

(Time taken will obviously vary between computers)

```
====================================
Degree 10
Slow Multiplication:   0.000036 secs
Fast Multiplication:   0.000465 secs
====================================
Degree 100
Slow Multiplication:   0.002695 secs
Fast Multiplication:   0.003318 secs
====================================
Degree 1000
Slow Multiplication:   0.274458 secs
Fast Multiplication:   0.032793 secs
====================================
Degree 10000
Slow Multiplication:  33.857816 secs
Fast Multiplication:   0.841920 secs
====================================
```
