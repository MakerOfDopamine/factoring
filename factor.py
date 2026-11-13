import math, random

def is_prime(n, k = 20):
    # Returns whether a number is prime.
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError("Invalid number!")

    if n <= 1: return False
    if n <= 3: return True

    if not (n & 1): return False

    # We use miller-rabin primality testing.
    # TODO actually implement this

    s = 0
    d = n - 1
    while not (d & 1):
        d >>= 1
        s += 1

    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        for _ in range(s):
            y = pow(x, 2, n)
            if (y == 1 and (x != 1 and x != (n - 1))):
                return False
            x = y
        if y != 1:
            return False

    return True

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

def pollard_rho(n, iters_max = 10000, k = 1):
    def g(x):
        return (pow(x, 2, n) + k) % n

    iters = 0
    a = 2
    b = a
    d = 1
    for _ in range(iters_max):
        a = g(a)
        b = g(g(b))
        d = math.gcd(abs(a - b), n)
        iters += 1

        if d != 1: 
            print(f"broke at iter {iters}")
            break

    if iters == iters_max:
        return None
    if d == n:
        print(f"found n")
        return pollard_rho(n, iters_max - iters, random.randint(1, n - 1))

    return d

def pollard_pm1(n, b1 = 1000, b2 = 0):
    '''Pollard's P-1 algorithm with Crandall's stage 2. Set b2 to 0 for auto, b2 to -1 to disable altogether.'''
    if not (n & 1): return 2

    M = sieve(b1)
    prods = []
    pp = 1
    for i, prime in enumerate(M):
        multiplicity = math.floor(math.log(b1) / math.log(prime))
        pp *= pow(prime, multiplicity)
        if (i + 1) % 1000 == 0:
            prods.append(pp)
            pp = 1
    prods.append(pp)
    #print(f"Calculated prod with {math.floor(math.log10(prod)) + 1} digits")
    g = 2
    for p in prods:
        g = pow(g, p, n)
        gcd_step = math.gcd(g + n - 1, n)
        if gcd_step != 1:
            if gcd_step == n:
                print("Found all factors of n! reduce b1 size")
                return None
            return gcd_step

    return None

def factor(n):
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError("Invalid number!")

    if n <= 1: return []
    elif n <= 3: return [n]

    # Pick factoring algorithms

def main():
    # Placeholder function. Eventually this wants the driver
    print(pollard_pm1(747867832862997038641, 150))

# boilerplate
if __name__ == "__main__":
    main()