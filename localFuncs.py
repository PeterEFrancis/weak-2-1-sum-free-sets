import numpy as np
# from addcomb import *
import itertools as it
import copy
import time
from IPython.display import display, Markdown, clear_output
from numba import njit


@njit
def D(n):
    output = np.array([n])
    for i in range(1,n):
        if (n % i) == 0:
            output = np.concatenate((output,np.array([i])))
    return output

@njit
def isPrime(n):
    for i in range(2,int(n**0.5)+1):
        if n%i==0:
            return False
    return True

@njit
def primeD(n):
    return np.array([i for i in D(n) if isPrime(i) and i!=1])


@njit
def has_pd23(n):
    for d in primeD(n):
        if d % 3 == 2:
            return True
    return False



def sumSet(A,B,mod):
    return {(a + b) % mod for a in A for b in B}

def check(A,B,C,mod):
    return len(sumSet(A,B,mod).intersection(C)) != 0

def sumSet2(A,mod):
    return {(a + b) % mod for a in A for b in A if a != b}

def check2(A,C,mod):
    return len(sumSet2(A,mod).intersection(C)) != 0

def stringSet(A):
    return ''.join(str(i) + ', ' for i in A)


def sumSet2R(A, zTup):
    sumSet = set()
    for a1 in A:
        for a2 in A:
            if ((a1[0] - a2[0]) % zTup[0] != 0) or ((a1[1] - a2[1]) % zTup[1] != 0):
                sumSet.add(((a1[0]+a2[0]) % zTup[0], (a1[1]+a2[1]) % zTup[1]))
    return sumSet


# looks for (2,1)-sum-free sets in Z_7 x Z_w
def look7(z:int, s:tuple, diff:int, verbose=False, enforceSameStart=False, number=-1):
    done = 1
    for a in range(z):
        A0 = {(a + diff * i) % z for i in range(s[0])}        
        if check2(A0,A0,z):
            continue
        for b in range(z):
            A1 = {(b + diff * i) % z for i in range(s[1])}
            if check(A0,A1,A1,z):
                continue
            for c in range(z):
                A6 = {(c + diff * i) % z for i in range(s[2])}
                if check(A0,A6,A6,z) or check(A1,A6,A0,z):
                    continue
                for d in range(z):
                    A2 = {(d + diff * i) % z for i in range(s[3])}
                    if check(A0,A2,A2,z) or check(A6,A2,A1,z) or check2(A1,A2,z):
                        continue
                    for e in range(z):
                        A5 = {(e + diff * i) % z for i in range(s[4])}
                        if check(A0,A5,A5,z) or check(A1,A5,A6,z) or check(A2,A5,A0,z) or check2(A6,A5,z):
                            continue
                        for f in range(z):
                            A3 = {(f + diff * i) % z for i in range(s[5])}
                            if check(A0,A3,A3,z) or check(A5,A3,A1,z) or check(A2,A3,A5,z) or check(A6,A3,A2,z):
                                continue
                            if check(A1,A2,A3,z) or check2(A5,A3,z) or check2(A3,A6,z):
                                continue
                            for g in range(z):
                                A4 = {(g + diff * i) % z for i in range(s[6])}
                                if  check(A0,A4,A4,z) or check(A1,A4,A5,z) or check(A3,A4,A0,z) or check(A5,A4,A2,z) or check(A2,A4,A6,z):
                                    continue
                                if check(A6,A4,A3,z) or check(A1,A3,A4,z) or check2(A4,A1,z) or check2(A2,A4,z) or check(A5,A6,A4,z):
                                    continue
                                    
                                # at this point, the set is good
                                
                                if not (enforceSameStart and not(a==b==c==d==e==f==g)):
    
                                    print(f"( {done} ) ----------------------------------------------------------------------" + ("*******" if (a==b==c==d==e==f==g) else ""))
                                    
                                    print(f"         a={a}, b={b}, c={c}, d={d}, e={e}, f={f}, g={g}        |A|={sum(len(A) for A in [A0,A1,A2,A3,A4,A5,A6])}\n")

                                    if verbose == True:
                                        A = [A0,A1,A2,A3,A4,A5,A6]
                                        if any(len(A[i]) != s[i] for i in range(7)):
                                            print("diff does not generate an A set of proper size\n")
                                            return False
                                        print(f" A0={mod_sort(z,A0)},\n A1={mod_sort(z,A1)},\n A2={mod_sort(z,A2)},\n A3={mod_sort(z,A3)},\n A4={mod_sort(z,A4)},\n A5={mod_sort(z,A5)},\n A6={mod_sort(z,A6)}\n")
                                        print('A = { ' + ''.join(''.join(f'({i},{j}), ' for j in [A0,A1,A2,A3,A4,A5,A6][i]) for i in range(7))[:-2] + ' }\n')

                                
                                if (number==-1) or done < number:
                                     done += 1
                                else:
                                    return True

    

    
@njit
def prod(tup):
    p = 1
    for el in tup:
        p *= el
    return p

                                
def mu21(G, verbose = False):
    if verbose:
        return mu(G,2,1,verbose = True)
    else:
        return int(v(1,G[-1],3)*prod(G)/(G[-1]))
                                
        
        
def isWeak21SumFree(A, zTup, verbose = False):
    if verbose:
        print(sumSet2R(A,zTup))
        print(" ")
        print(sumSet2R(A,zTup).intersection(A))
    return sumSet2R(A,zTup).intersection(A) == set()







