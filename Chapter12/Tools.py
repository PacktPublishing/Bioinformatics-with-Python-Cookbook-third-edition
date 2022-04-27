import functools


def fibo_iter(n):
    if n == 0:
        return 0
    if n == 1:
        return 1
    last = 1
    second_last = 1
    for i in range(3, n + 1):
        result = second_last + last
        second_last = last
        last = result
    return result


def fibo_naive(n):
    if n == 0:
        return 0
    if n == 1:
        return 1
    return fibo_naive(n - 1) + fibo_naive(n - 2)


@functools.lru_cache
def fibo(n):
    if n == 0:
        return 0
    if n == 1:
        return 1
    return fibo(n - 1) + fibo(n - 2)


time fibo_iter(100)
#time fibo_naive(1000)
time fibo(1000)


def factorial(n):
    if n == 1:
        return 1
    return n * factorial(n - 1)


factorial(20000)
