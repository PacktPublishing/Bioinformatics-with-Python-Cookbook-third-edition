def fibo_iter(n):
    if n < 2:
        return n
    last = 1
    second_last = 0
    for _i in range(2, n + 1):
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


fibo_iter(0)
fibo_iter(1)
fibo_iter(2)
fibo_iter(3)
fibo_iter(4)
fibo_iter(5)
fibo_iter(6)
fibo_naive(1000)


def factorial(n):
    if n == 1:
        return 1
    return n * factorial(n - 1)


factorial(5)
factorial(20000)
