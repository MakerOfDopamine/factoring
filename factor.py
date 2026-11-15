import math, random, shutil, os, time

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

def fermat(n, rounds = 100_000):
    a = math.isqrt(n)
    if a * a == n: return a
    for _ in range(rounds):
        a += 1
        b = a * a - n
        bs = math.isqrt(b)
        if bs * bs == b:
            return a - bs

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
            #print(f"broke at iter {iters}")
            break

    if iters == iters_max:
        return None
    if d == n:
        #print(f"found n")
        return pollard_rho(n, iters_max - iters, random.randint(1, n - 1))

    return d

def pollard_pm1(n, b1 = 1000, b2 = 0):
    '''Pollard's P-1 algorithm with Crandall's stage 2. Set b2 to 0 for auto, b2 to -1 to disable altogether.'''
    if b2 == 0: b2 = b1 * 100
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
                #print("Found all factors of n! reduce b1 size")
                return None
            return gcd_step

    # part 2
    if b2 == -1: return None
    M2 = seg_sieve(b1, b2)
    if len(M2) == 0: return None
    largest_diff = max([M2[i] - M2[i-1] for i in range(1, len(M2) - 1)])
    powers = {}
    base = pow(g, 2, n)
    last = 1
    for i in range(largest_diff // 2):
        last *= base
        last %= n
        powers[(i + 1) * 2] = last

    z = pow(g, M2[0], n)
    acc = (z - 1) % n
    for p, q in zip(M2, M2[1:]):
        z = z * powers[q - p] % n
        acc = acc * (z - 1) % n
     
    gcd_step = math.gcd(acc, n)
    if gcd_step != 1:
        if gcd_step == n:
            #print("Found all factors of n! reduce b1 size")
            return None
        return gcd_step
    return None

def ecm(n, b1 = 1000, b2 = 0, curves = 1):
    if n % 2 == 0: return 2
    if n % 3 == 0: return 3
    if n % 5 == 0: return 5

    primes = sieve(b1)
    primes2 = seg_sieve(b1, b2)

    #D = 210 if b2 < 1000000 else 2310
    D = 30
    # suyama's parameterization
    for curve in range(curves):
        #print("curve", curve)
        sigma = random.randint(6, n - 2)
        u = (sigma * sigma - 5) % n
        v = (sigma * 4) % n

        x0 = pow(u, 3, n)
        z0 = pow(v, 3, n)

        try:
            # because the inverse can fail for some reason
            a = (pow(v - u, 3, n) * (3 * u + v) * pow(4 * u * u * u * v, -1, n) - 2) % n
            a24 = ((a + 2) * pow(4, -1, n)) % n
        except ValueError:
            factor = math.gcd(4 * u * u * u * v, n)
            if factor == n:
                continue
            else:
                return factor
        
        def _elliptic_add(x0, x1, z0, z1, d0, d1):
            # Why must elliptic curves have insane rules?
            U = ((x0 - z0) * (x1 + z1)) % n
            V = ((x0 + z0) * (x1 - z1)) % n

            return (d1 * pow(U + V, 2, n)) % n, (d0 * pow(U - V, 2, n)) % n
        
        def _elliptic_double(x, z):
            s = pow(x + z, 2, n)
            d = pow(x - z, 2, n)
            t = s - d
            return (s * d) % n, (t * (d + a24 * t)) % n

        def _elliptic_mul(Px, Pz, k):
            pax, paz = Px, Pz
            pbx, pbz = _elliptic_double(Px, Pz)

            if k == 1: return pax, paz

            for bit in bin(k)[3:]:
                if bit == "0":
                    pbx, pbz = _elliptic_add(pax, pbx, paz, pbz, Px, Pz)
                    pax, paz = _elliptic_double(pax, paz)
                else:
                    pax, paz = _elliptic_add(pax, pbx, paz, pbz, Px, Pz)
                    pbx, pbz = _elliptic_double(pbx, pbz)

            return pax, paz


        Qx, Qz = x0, z0
        for p in primes:
            exp = math.floor(math.log(n) / math.log(p))
            Qx, Qz = _elliptic_mul(Qx, Qz, pow(p, exp))

        g = math.gcd(Qz, n)
        if g != 1:
            if g == n:
                pass
                #print("Found n")
            else:
                return g

        # part 2

        j = {}
        for i in range(1, D // 2, 2):
            if math.gcd(i, D) == 1:
                j[i] = tuple(_elliptic_mul(Qx, Qz, i))

        Sx, Sz = _elliptic_mul(Qx, Qz, D)
        m = max(math.floor(b1 / D), 2)
        Rx, Rz = _elliptic_mul(Sx, Sz, m)
        Rpx, Rpz = _elliptic_mul(Sx, Sz, m - 1)

        acc = 1
        #print("n sigma j m", n, sigma, j, m)
        #print("S", Sx, Sz)
        #print("R", Rx, Rz)
        #print("Rp", Rpx, Rpz)
        while (m * D - D//2) <= b2:
            for jval in j:
                if (is_prime(m * D + jval) and (b1 < m * D + jval <= b2)) or (is_prime(m * D - jval) and (b1 < m * D - jval <= b2)):
                    acc = (acc * (Rx * j[jval][1] - Rz * j[jval][0])) % n
                    #print("m", m, "j", jval, "acc", acc)
            Rx, Rz, Rpx, Rpz = *_elliptic_add(Rx, Sx, Rz, Sz, Rpx, Rpz), Rx, Rz

            g = math.gcd(acc, n)
            if g != 1:
                if g == n:
                    pass
                else: return g
            m += 1

    return None

def _sym_legendre(n, p):
    # Returns the legendre symbol (N : p).
    if not is_prime(p):
        print("p is not prime")
        return None
    n %= p
    if n == 0: return 0
    a = pow(n, (p - 1) // 2, p)
    if a == 1:
        return a
    else:
        return -1

def sqrt_mod(n, p):
    # Returns n^(1/2) mod p.
    if not is_prime(p):
        print("p is not prime")
        return None
    if _sym_legendre(n, p) == -1:
        print(n, "has no sqrt", p)
        return None
    if p == 2:
        return n % 1

    n %= p
    q = p - 1
    s = 0
    while not (q & 1):
        q //= 2
        s += 1
    z = 3
    while _sym_legendre(z, p) != -1:
        z += 1
    M = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    R = pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return R
        tt = t
        i = 0
        while tt != 1: 
            tt = pow(tt, 2, p)
            i += 1
        b = pow(c, 1 << (M - i - 1), p)
        M = i
        c = pow(b, 2, p)
        t = (c * t) % p
        R = (R * b) % p

def qs(n):
    # Quadratic Sieve.
    B = 1000
    primes = sieve(B)
    sieve_primes = []
    for p in primes:
        if _sym_legendre(n, p) == 1:
            sieve_primes.append(p)
    

def factor(n):
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError("Invalid number!")

    if n <= 1: return []
    elif n <= 3: return [n]

    # Pick factoring algorithms
    factors = []
    f = trial_factor(n, min(math.isqrt(n), 1_000_000))
    factors.extend(f[0])
    n = f[1]

    if (is_prime(n)):
        factors.append(n)
        n = 1
        
    if n == 1:
        #print(f"{math.prod(factors)} = {" x ".join(map(str, factors))}")
        return (factors, n)

    fa = 0
    while fa is not None:
        fa = pollard_rho(n, 1000000)
        if fa is not None:
            factors.append(fa)
            n //= fa
        if n == 1:
            #print(f"{math.prod(factors)} = {" x ".join(map(str, factors))}")
            return (factors, n)

    fa = pollard_pm1(n, 750, 75000)
    if fa is not None:
        factors.append(fa)
        n //= fa
    if n == 1:
        return (factors, n)

    for i in range(25):
        fa = ecm(n, 200, 20000)
        if fa is not None:
            factors.append(fa)
            n //= fa
        if n == 1:
            return (factors, n)

    fa = pollard_pm1(n, 5000, 500000)
    if fa is not None:
        factors.append(fa)
        n //= fa
    if n == 1:
        return (factors, n)

    for i in range(50):
        fa = ecm(n, 500, 50000)
        if fa is not None:
            factors.append(fa)
            n //= fa
        if n == 1:
            return (factors, n)

    return (factors, n)

def main():
    # Placeholder function. Eventually this wants the driver
    #print(pollard_pm1(191919191919191919191919191919191917 // 38567, 1000000, 100000000))
    #print(ecm(270096654252030534964291310083834119193, 100000, 0, 100))
    #exit()
    os.system("cls" if os.name == "nt" else "clear")
    BASIC = """
=====================================================
| Welcome to Factr!                                 |
| (1) Auto factor input number                      |
| (2) Run Trial Factoring on input number           |
| (3) Run Fermat's Factorization Method             |
| (4) Run Pollard's Rho                             |
| (5) Run Pollard's P-1                             |
| (6) Run William's P+1                             |
| (7) Run Elliptic Curve Method                     |
| (8) Run Quadratic Sieve                           |
| (9) Open miscellaneous factorization methods menu |
| (0) EXIT                                          |
=====================================================
        """
    while True:
        print("\033[H", end="", flush=True)
        print(BASIC)
        option = input("> ")
        match option:
            case "1":
                N = int(input("Input your number > "))
                f, cf = factor(N)
                print(f"{math.prod(f) * cf} = {" x ".join(map(str, f))}{f" x {"C!" if not is_prime(cf) else ""}{cf}" if cf > 1 else ""}")
            case "2":
                N = int(input("Input your number > "))
                x = int(input("Input trial factor bound > "))
                f, cf = trial_factor(N, x)
                if len(f) == 0:
                    print("Found no factors.")
                else:
                    print(f"{N} = {" x ".join(map(str, f))}{f" x {"C!" if not is_prime(cf) else ""}{cf}" if cf > 1 else 0}")
            case "0":
                print("Thanks for using Factr!")
                return
            case _:
                print("That's not a valid option!")
        input("continue... ")

        print("\033[H", end="", flush=True)
        print(f"{" " * shutil.get_terminal_size().columns}\n" * (shutil.get_terminal_size().lines - 1), end="")
        print("\033[H", end="", flush=True)

# benchmarking suite
def bench_pm1_b1(ss = 128, r = 32):
    # I originally planned on doing this but realized it took too long and was messy
    for size in range(ss, ss + r):
        print("testing size", size, "bit numbers")
        found = [0] * 20
        times = [0] * 20
        num = 0
        while num < 100:
            print(num, end="\r")
            n = trial_factor(random.randint(1 << size, 1 << (size + 1)), 100_000)[1]
            if n < (1 << size): continue
            least_findable = 0
            for steps in range(1, 21):
                t = time.time()
                if pollard_pm1(n, 50 * steps, -1): 
                    least_findable = steps - 1
                    times[steps - 1] += time.time() - t
                    break
                times[steps - 1] += time.time() - t

            for i in range(least_findable, 20):
                found[i] += 1
            num += 1

        rate = [i / j for i, j in zip(found, times)]
        opt = rate.index(max(rate))
        print("Optimal for bit", size, "is B1", 50 * (opt + 1), "found", found[opt], "times")

def bench_pm1_b1_v2(ss = 128, r = 32):
    for size in range(ss, ss + r):
        #print("Testing", size, "bit numbers")
        b1s = []
        b2s = []
        num = 0
        while num < 200:
            #print(num, end="\r")
            p = 4
            q = 4
            while not is_prime(p):
                p = random.randint(1 << (size // 2), 2 << (size // 2))
            while not is_prime(q):
                q = random.randint(1 << ((size + 1) // 2), 2 << ((size + 1) // 2))

            if p > q: p, q = q, p
            n = p * q
            pm1 = p - 1
            factors, core = trial_factor(pm1, 100_000)
            if is_prime(core):
                factors.append(core)
                core = 1
            while core > 1:
                t = ecm(core, 100, 0, 20)
                if t is not None:
                    core //= t
                    factors.append(t)
                else: break

            if core > 1: continue
            factors.sort()
            req_b2 = factors[-1] + 1
            req_b1 = factors[-2] + 1

            if req_b2 > req_b1 * 100: continue

            b1s.append(req_b1)
            b2s.append(req_b2)
            num += 1

        #print("Required B1:", round(sum(b1s) / len(b1s)), "Required B2:", round(sum(b2s) / len(b2s)))
        #print(size, round(sum(b1s) / len(b1s)), round(sum(b2s) / len(b2s)))

        print("\n".join(map(str, b1s)))
        print()
        print("\n".join(map(str, b2s)))

def bench_pm1_b2():
    for size in range(48, 64):
        #print("testing size", size, "bit numbers")
        found = [0] * 20
        num = 0
        while num < 100:
            #print(num, end="\r")
            n = trial_factor(random.randint(1 << size, 1 << (size + 1)), 1_000_000)[1]
            if n < (1 << (size - 4)): continue
            least_findable = 0
            for steps in range(1, 21):
                if pollard_pm1(n, 500, 500 * steps): 
                    least_findable = steps - 1
                    break

            for i in range(least_findable, 20):
                found[i] += 1
            num += 1
        print(found)

def bench2():
    # benchmarking function to find crossovers between p-1 and ecm.
    pm1_s1 = 0
    pm1_s2 = 0
    pm1_s3 = 0

    for i in range(100):
        print("Testing number", i+1)
        n = trial_factor(random.randint(2 << 64, 2 << 65), 1000000)[1]
        if pollard_pm1(n, 1000):
            pm1_s1 += 1
            pm1_s2 += 1
            pm1_s3 += 1
        elif pollard_pm1(n, 10000):
            1

# boilerplate
if __name__ == "__main__":
    #main()
    #bench_pm1_b1_v2(96, 1)
    #print("factor", ecm(183711993843743, 50, 250))
    print(fermat(51983))