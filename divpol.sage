def division_poly(l, x, A, B):
    f = x^3 + A*x + B
    memo = {}
    
    def g(k):
        if k in memo:
            return memo[k]
        if k == 0: res = x.parent()(0)
        elif k in (1, 2): res = x.parent()(1)
        elif k == 3: res = 3*x^4 + 6*A*x^2 + 12*B*x - A^2
        elif k == 4: res = 4*x^6 + 20*A*x^4 + 80*B*x^3 - 20*A^2*x^2 - 16*A*B*x - 4*A^3 - 32*B^2
        else:
            m = k // 2
            if k % 2 == 1:  # Odd k
                if m % 2 == 1:
                    res = g(m+2) * g(m)^3 - 16 * f^2 * g(m-1) * g(m+1)^3
                else:
                    res = 16 * f^2 * g(m+2) * g(m)^3 - g(m-1) * g(m+1)^3
            else:  # Even k
                res = g(m) * (g(m+2) * g(m-1)^2 - g(m-2) * g(m+1)^2)
        memo[k] = res
        if k==2:
            return f
        else:
            return res

    return [g(i) for i in l]


#primes = primes_first_n(10)
#primes = primes_first_n(6)
#primes.extend([37])
primes=[2,3,5]
n = len(primes)
p = product(primes) - 1
is_prime(p)
F=GF(p)
R = PolynomialRing(F, 2*n + 1,'X')
X=R.gens()
L=division_poly(primes,x=X[0],A=X[2*n], B=-1-X[2*n])
for i in range(n):
    L[i] = L[i].subs({X[0]:X[i]})
    L.extend([X[n+i]^2-X[i]^3-X[2*n]*X[i]+1+X[2*n]])
I=R.ideal(L)
I
I.groebner_basis()
decomp = I.primary_decomposition()
