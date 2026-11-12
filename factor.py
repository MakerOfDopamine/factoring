import math

def is_prime(n, k = 20):
    # Returns whether a number is prime.
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError("Invalid number!")

    if n <= 1: return False
    if n <= 3: return True

    # We use miller-rabin primality testing.
    # TODO actually implement this

def sieve(n):
    # Returns primes up to n.
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError("Invalid number!")

    if n <= 1: return []
    if n == 2: return [2]
    if n == 3: return [2, 3]

    primes = bytearray([1]) * (n + 1)
    primes[0] = 0
    primes[1] = 0
    pt = 2
    while pt < math.isqrt(n) + 1:
        if primes[pt]:
            primes[pt * pt :: pt] = bytes(len(range(pt * pt, n+1, pt)))
        pt += 1
    return [val for val, isp in enumerate(primes) if isp]

def seg_sieve(n1, n2):
    # Returns primes between n1 and n2 (range conventions, including n1 but not including n2).
    if not isinstance(n1, int) or isinstance(n1, bool):
        raise TypeError("Invalid number!")
    if not isinstance(n2, int) or isinstance(n2, bool):
        raise TypeError("Invalid number!")

    primes = bytearray([1]) * (n2 - n1)
    for p in sieve(math.isqrt(n2)):
        start = max(p * p, (n1 + p-1) // p * p)
        primes[(start - n1) :: p] = bytes(len(range(start, n2, p)))
    return [val + n1 for val, isp in enumerate(primes) if isp]

# the 1_000_000 notation used in the following code is convenience because I do not trust my ability to Type Zeroes Correctly
def trial_factor(n, lim = 0):
    factors = []
    if (lim == 0):
        lim = math.isqrt(n)

    # if lim is small enough we do not need segmented sieve.
    if (lim < 10_000_000):
        primes = sieve(lim)
        for p in primes:
            while n % p == 0:
                n //= p
                factors.append(p)
                if n == 1:
                    break

    # segmented sieve 
    elif (lim < 10_000_000_000):
        for i in range(0, lim, 10_000_000):
            print(i)
            p_div = seg_sieve(max(i, 2), min(i + 10_000_000, lim))
            for p in p_div:
                while n % p == 0:
                    n //= p
                    factors.append(p)
                    if n == 1:
                        break

    # sanity check so lim isn't set too high
    else:
        print("Trial factor threshold too high! Consider using pollard's rho or ECM.")

    return (factors, n)

def factor(n):
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError("Invalid number!")

    if n <= 1: return []
    elif n <= 3: return [n]

    # Pick

def main():
    # Placeholder function. Eventually this wants the driver
    print(trial_factor(260007297762127))

# boilerplate
if __name__ == "__main__":
    main()