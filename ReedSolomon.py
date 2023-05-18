from gmpy2 import mpz, next_prime
from math import ceil
import random

def binaryEGCD(a, b):
    r = mpz(a)
    rr = mpz(b)
    e = 0
    while (r % 2 == 0 and rr % 2 == 0):
        r >>= 1
        rr >>= 1
        e += 1
    
    aa = r
    bb = rr
    s = mpz(1)
    t = mpz(0)
    ss = mpz(0)
    tt = mpz(1)

    while (rr != 0):
        while (r % 2 == 0):
            r >>= 1
            if (s % 2 == 0 and t % 2 == 0):
                s >>= 1
                t >>= 1
            else:
                s = (s + bb) >> 1
                t = (t - aa) >> 1
            
        while (rr % 2 == 0):
            rr >>= 1
            if (ss % 2 == 0 and tt % 2 == 0):
                ss >>= 1
                tt >>= 1
            else:
                ss = (ss + bb) >> 1
                tt = (tt - aa) >> 1
        
        if (rr < r):
            (r, s, t, rr, ss, tt) = (rr, ss, tt, r, s, t)
        rr -= r
        ss -= s
        tt -= t
    
    d = (1 << e) * r
    return (d, s, t)

def EGCD(a, b):
    r = [mpz(a), mpz(b)]
    s = [mpz(1), mpz(0)]
    t = [mpz(0), mpz(1)]

    while (1):
        qi = a // b # these are previous ones
        ri = a % b
        si = s[-2] - qi * s[-1] # s(i+1) = s(i-1) - q(i)s(i)
        ti = t[-2] - qi * t[-1] # t(i+1) = t(i-1) - q(i)t(i) 

        r.append(mpz(ri))
        s.append(mpz(si))
        t.append(mpz(ti))
        a = b
        b = ri
        if (ri == 0):
            break

    return (r, s, t)

def inverse(b, n): #b^(-1) mod(n)
    (d, s, t) = binaryEGCD(n, b)
    if (d != 1): # inverse doesn't exist
        return -1
    elif (t < 0):
        return n + t
    else:
        return t

def CRT(A, N):
    n = mpz(1)
    for x in N:
        n *= mpz(x)
    a = 0

    for i in range (0, len(A)):
        ni = mpz(n) // mpz(N[i])

        # brute force method to calculate ni
        # ni = mpz(1)
        # for j in range (0, len(N)):
        #     if (i == j):
        #         continue
        #     ni *= mpz(N[j])

        bi = ni % N[i]
        ti = inverse(bi, N[i])
        ei = (mpz(ni) * mpz(ti)) % n
        a = (a + (mpz(A[i]) * ei) % n) % n

    return a

class ReedSolomon:
    def __init__(self, M, mu): # functions as GlobalSetup(mu, M)
        self.M = M
        self.mu = mu
        self.N = [2] # ni's
        mulPrefix = [1, mpz(2)] # mulPrefix[i] stores product of N[j] such that j < i 
        self.P = mpz(2) # since k is atleast 1
        self.n = mpz(2) # since k is atleast 1
        while (self.n <= 2 * self.P * self.P * self.M):
            self.N.append(next_prime(self.N[-1]))
            self.n *= self.N[-1]
            mulPrefix.append(self.n)
            l = mpz(ceil(mu * len(self.N)))
            self.P = mulPrefix[-1] // mulPrefix[len(mulPrefix) - l - 1]
            if (self.N[-1] > (1 << 16)): # bound on the max prime number taken - equivalent to bound on k
                print(len(self.N))
                break

            # bruteforce method to calculate P
            # self.P = mpz(1)
            # for i in range(l):
            #     self.P *= self.N[-1 - i]

        self.k = len(self.N)

    def ReedSolomonSend(self, a):
        A = []
        for x in self.N:
            A.append(a % x)
        return self.Transmit(A)
    
    def Transmit(self, A):
        l = random.randint(0, ceil(self.mu * self.k)) # random number of positions which will get corrupted
        map = []
        for i in range(self.k):
            map.append(0)
        ll = 0
        B = []
        while (ll < l):
            rn = random.randint(0, self.k - 1) # randomly choosing which l positions will get corrupted
            if (map[rn] == 0):
                map[rn] = 1
                ll += 1
        for i in range(0, self.k):
            if (map[i] == 1):
                bi = random.randint(0, self.N[i] - 1) # corrupting 
                while (bi == A[i]):
                    bi = random.randint(0, self.N[i] - 1)
                B.append(bi)
            else:
                B.append(A[i])
        
        return B
    
    def ReedSolomonReceive(self, B):
        b = CRT(B, self.N) 
        rr = mpz(self.M * self.P) # r* = MP in rational reconstruction
        ri, si, ti = EGCD(self.n, b)
        for j in range(0, len(ri)):
            if (ri[j] <= rr): 
                rd = ri[j]
                td = ti[j]
                break

        if (td == 0 or rd % td != 0):
            return -1
        else:
            a = (rd // td)
            return a
        
def nonNegativeMu():
    mu = float(input("Enter corruption factor: "))
    while (mu < 0):
        print("Corruption factor can't be negative. Retry")
        mu = float(input("Enter corruption factor: "))
    return mu

def inputvalues():
    M = int(input("Enter M: "))
    while (M < 0):
        print("M can't be negative. Retry")
        M = int(input("Enter M: "))
    mu = nonNegativeMu() 
    while (mu > 0.46):
        print("Warning! The corruption factor is too large and it might not be able reconstruct messages.")
        print("Press 0: Continue")
        print("Press 1: Enter new corruption factor")
        choice = int(input())
        if (choice == 0):
            break
        elif (choice == 1):
            mu = nonNegativeMu()
        else:
            print("Enter valid choice.")

    return (M, mu)

M, mu = inputvalues()
rs = ReedSolomon(M, mu)
while (1):
    print("Press 1: To enter new message to transmit")
    print("Press 2: Change (M, mu)")
    print("Press 3: Exit")
    choice = int(input())
    if (choice == 1):
        a = int(input("Enter message bounded by M: "))
        while (a > M):
            print("Message must be less than or equal to M.")
            a = int(input("Enter message bounded by M: "))
            
        B = rs.ReedSolomonSend(a)
        b = rs.ReedSolomonReceive(B)
        print(f'Reconstructed message: {b}')
        if (a == b):
            print("Successfully reconstructed.")
        else:
            print("Reconstruction unsuccessful.")
    elif (choice == 2):
        M, mu = inputvalues()
        rs = ReedSolomon(M, mu)
    elif (choice == 3):
        break
    else:
        print("Enter valid choice.")

