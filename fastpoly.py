#!/usr/bin/env python3

import math
import cmath
import random
import time

class Polynomial:
    def __init__(self, coeff=[0]):
        while coeff[-1] == 0 and len(coeff) > 1:
            coeff.pop()
        self.coeff = coeff

    def __repr__(self):
        output = "{0:d}".format(self.coeff[0])
        if len(self.coeff) > 1:
            term = self.coeff[1]
            if term > 0:
                output += " + "
            elif term < 0:
                output += " - "
                term *= -1
            if term != 0:
                output += "{0:d} x".format(term)
        for exponent in range(2, len(self.coeff)):
            term = self.coeff[exponent]
            if term > 0:
                output += " + "
            elif term < 0:
                output += " - "
                term *= -1
            else:
                continue

            output += "{0:d} x^{1}".format(term, exponent)

        return output

def slowPolyMult(p, q):
    degp = len(p.coeff) - 1
    degq = len(q.coeff) - 1
    pqcoeff = [0] * (degp + degq + 1)
    for d in range(degp + degq + 1):
        for i in range(max(0, d-degq), min(d+1, degp + 1)):
            pqcoeff[d] += p.coeff[i] * q.coeff[d-i]
    return Polynomial(pqcoeff)

def FFT(p, omegas):
    # Convert the p as a list of coefficients                 #
    # to the value representation using FFT with nth roots of #
    # unity in omegas where n > deg p                         #
    n = len(omegas)
    if n == 1:
        return p
    even = FFT(p[::2], omegas[::2])
    odd = FFT(p[1::2], omegas[::2])
    values = [0] * n
    for i in range(n // 2):
        values[i] = even[i] + omegas[i] * odd[i]
        values[i + n // 2] = even[i] - omegas[i] * odd[i]
    return values

def IFFT(p, omegas):
    # Convert the p as a list of values to the coefficient #
    # representation using FFT with nth roots of unity in  #
    # omegas where n > deg p                               #
    n = len(omegas)
    if n == 1:
        return p
    even = IFFT(p[::2], omegas[::2])
    odd = IFFT(p[1::2], omegas[::2])
    values = [0] * n
    for i in range(n // 2):
        values[i] = 0.5 * (even[i] + omegas[i] * odd[i])
        values[i + n // 2] = 0.5 * (even[i] - omegas[i] * odd[i])
    return values

def fastPolyMult(p, q):
    pco = p.coeff
    qco = q.coeff
    degp = len(pco) - 1
    degq = len(qco) - 1
    fftdeg = 2**math.ceil(math.log(degp + degq + 1, 2))
    omegas = [0] * fftdeg
    invomegas = [0] * fftdeg
    for i in range(fftdeg):
        omegas[i] = cmath.exp(i * cmath.pi * 2j / fftdeg)
        invomegas[i] = cmath.exp(i * cmath.pi * -2j / fftdeg)
    pco += [0] * (fftdeg - degp - 1)
    qco += [0] * (fftdeg - degq - 1)
    pvalues = FFT(pco, omegas)
    qvalues = FFT(qco, omegas)
    pqvalues = [0]*fftdeg
    for i in range(fftdeg):
        pqvalues[i] = pvalues[i] * qvalues[i]
    pqcoeff = [round(c.real) for c in IFFT(pqvalues, invomegas)]
    return Polynomial(pqcoeff)

def randomPoly(degree):
    coeff = [0] * (degree + 1)
    for i in range(degree + 1):
        coeff[i] = random.randint(-10, 11)
    return Polynomial(coeff)



ROUNDS = 10
slowtime = 0
fasttime = 0

print("=" * 36)

for deg in [10,100,1000,10000]:
    for i in range(ROUNDS):
        p = randomPoly(deg)
        q = randomPoly(deg)

        start = time.perf_counter()
        a = slowPolyMult(p, q)
        end = time.perf_counter()
        slowtime += (end - start)

        start = time.perf_counter()
        b = fastPolyMult(p, q)
        end = time.perf_counter()
        fasttime += (end-start)

        if (a.coeff != b.coeff):
            print("Error: Products not equal")
            break

    print("Degree {0:d}".format(deg))
    print("Slow Multiplication: {0:10.6f} secs".format(slowtime / ROUNDS))
    print("Fast Multiplication: {0:10.6f} secs".format(fasttime / ROUNDS))
    print("=" * 36)
