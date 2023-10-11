from constants import *

"""Discrete Fourier Transform for an array of values x."""
def dft(x: list) -> list:
    len_x = len(x)
    base = -2j * PI / len_x
    result = []
    for k in range(len_x):
        term = 0
        for n, x_n in enumerate(x):
            exp = base * k * n
            term += x_n * (EULER ** exp)
        result.append(term)

    return result

"""Fast Fourier Transform for an array of values x."""
def fft(x: list) -> list:
    len_x = len(x)
    if len_x % 2 > 0:
        raise ValueError("Array size must a Power of 2")
    elif len_x <= 2:
        return dft(x)
    else:
        base = -2j * PI / len_x
        evens = fft(x[::2])
        odds = fft(x[1::2])
        terms = [EULER ** (base * k) for k in range(len_x)]

        left = []
        right = []
        for i in range(int (len_x / 2)):
            l_term = evens[i] + terms[i] * odds[i]
            s = int(len_x / 2) + i
            r_term = evens[i] + terms[s] * odds[i]
            left.append(l_term)
            right.append(r_term)

        left.extend(right)
        return left
