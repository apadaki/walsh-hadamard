import numpy as np
import math, random
import pickle
from numpy.random.mtrand import randint
from scipy.linalg import lu, pinv

from computations import iterative_comp

# the standard flip factorization
def flip_factorization(H, k, soln):
    # for a solution L, express H_k = L + S
    L = soln

    S = H[k]-L

    # decompose L = uv^T
    v = np.empty((2**k, 1), dtype=int)
    u = np.empty((2**k, 1), dtype=int)
    v[:,0] = L[0]
    for i in range(2**k):
        u[i][0] = 1 if L[i][0] == L[0][0] else -1

    # factorize H_k as (A x B) using identity matrix concatenation
    A = np.concatenate((np.identity(2**k), u), axis=1)
    B = np.concatenate((S, v.T), axis=0)

    return A,B

# the identity factorization
def identity_factorization(H,k):
    return np.identity(2**k), H[k]

# a bunch of other factorizations
def flip_factorization_QR(H, k, soln):
    # for a solution L, express H_k = L + S
    L = soln

    S = H[k]-L

    # decompose L = uv^T
    v = np.empty((2**k, 1), dtype=int)
    u = np.empty((2**k, 1), dtype=int)
    v[:,0] = L[0]
    for i in range(2**k):
        u[i][0] = 1 if L[i][0] == L[0][0] else -1

    # UNIQUE STEP: factorize S = (QS x RS) with QR-decomposition
    QS, RS = np.linalg.qr(S)
    A = np.concatenate((QS, u), axis=1)
    B = np.concatenate((RS, v.T), axis=0)

    return A,B

def flip_factorization_LU(H, k, soln):
    # for a solution L, express H_k = L + S
    L = soln

    S = H[k]-L

    # decompose L = uv^T
    v = np.empty((2**k, 1), dtype=int)
    u = np.empty((2**k, 1), dtype=int)
    v[:,0] = L[0]
    for i in range(2**k):
        u[i][0] = 1 if L[i][0] == L[0][0] else -1

    # UNIQUE STEP: factorize S = (PxL) x Y with LU-decomposition
    P,L,Y = lu(S)
    A = np.concatenate((P@L, u), axis=1)
    B = np.concatenate((Y, v.T), axis=0)

    return A,B

def QR_factorization(H,k):
    Q,R = np.linalg.qr(H[k])
    return Q,R

def LU_factorization(H,k):
    p,l,u = lu(H[k])
    return p@l,u

def evaluate_factorization(A,B,k):
    # find number of nonzero elements in kronecker product of A,B^T
    nz = np.count_nonzero(A) * np.count_nonzero(B)
    return math.log2(nz)/(2*k)


def run_factorizations(H):
    for k in range(1,5):
        solutions_k = iterative_comp(H, k, verbose=False)

        sym_soln = solutions_k[0]
        for soln in solutions_k:
            if (soln!=soln.T).all():
                print('[non symmetric found]')
                sym_soln = soln
                break

        print('RESULTS (k = {}):'.format(k))
        A,B = identity_factorization(H,k)
        assert(np.allclose(A@B, H[k]))
        e_value = evaluate_factorization(A,B,k)
        print('identity factorization: ', e_value)

        A,B = flip_factorization(H, k, sym_soln)
        assert(np.allclose(A@B, H[k]))
        e_value = evaluate_factorization(A,B,k)
        print('flip factorization: ', e_value)

        A,B = flip_factorization_QR(H, k, sym_soln)
        assert(np.allclose(A@B, H[k]))
        e_value = evaluate_factorization(A,B,k)
        print('flip factorization (QR): ', e_value)

        A,B = flip_factorization_LU(H, k, sym_soln)
        assert(np.allclose(A@B, H[k]))
        e_value = evaluate_factorization(A,B,k)
        print('flip factorization (LU): ', e_value)

        A,B = QR_factorization(H,k)
        assert(np.allclose(A@B, H[k]))
        e_value = evaluate_factorization(A,B,k)
        print('QR factorization: ', e_value)

        A,B = LU_factorization(H,k)
        assert(np.allclose(A@B, H[k]))
        e_value = evaluate_factorization(A,B,k)
        print('LU factorization: ', e_value)

        print()

# generate a random sparse candidate for left factor and find e-value
def brute_force_bernoulli(H,k):
    A = np.matrix(np.zeros((2**k, 2**k+random.randint(1,3))))
    num_entries_A = A.shape[0]*A.shape[1]
    desired_nz_A = 2**k*0.9
    # choose Bernoulli parameter based on desired sparsitys
    p = desired_nz_A/num_entries_A
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            r = random.random()
            if r < p:
                A[i,j] = random.randint(0,1)*2-1
            else:
                A[i,j] = 0

    # find Moore-Penrose pseudo-inverse of A to compute B
    Ainv = pinv(A)

    B = np.round(Ainv@H[k], decimals=4)

    # verify that H[k] = AB, which is not necessarily given
    if np.allclose(A@B, H[k]):
        return A,B,evaluate_factorization(A,B,k)
    return -1,-1,-1


if __name__ == '__main__':
    # construct H1 through H7
    H = [np.array([1], dtype=int)]
    h1 = np.array([[1,1],[1,-1]])
    for i in range(1,8):
        H.append(np.kron(h1,H[i-1]))

    # run_factorizations(H)

    min_A,min_B=None,None
    min_e = 10
    n_iter = 100000
    for i in range(n_iter):
        if (i != 0 and i % 10000 == 0):
            print(i)
        A,B,e = brute_force_bernoulli(H,3)
        if e != -1 and e < min_e:
            min_A = A
            min_B = B
            min_e = e
    print(min_A)
    print(min_B)
    print(min_e)
