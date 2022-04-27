import functools


@functools.cache
def fibo(n):
    if n == 0:
        return 0
    if n == 1:
        return 1
    return fibo(n - 1) + fibo(n - 2)


fibo(1000)


def gene_min_reads(source, min_reads):
    return map(
        lambda x: x[0],
        filter(
            lambda x: x[1] >= min_reads,
            source.items()))


list(gene_min_reads({'LCT': 10, 'MRAP2': 1}, 2))


multiplication = lambda x, y: x * y

double = functools.partial(multiplication, 2)

double(3)
