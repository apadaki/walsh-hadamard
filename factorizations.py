import numpy as np
import math
import pickle
from scipy.linalg import lu

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

if __name__ == '__main__':
    # construct H1 through H7
    H = [np.array([1], dtype=int)]
    h1 = np.array([[1,1],[1,-1]])
    for i in range(1,8):
        H.append(np.kron(h1,H[i-1]))
    
    for k in range(1,5):
        solutions_k = iterative_comp(H, k, verbose=False)

        sym_soln = solutions_k[0]
        for soln in solutions_k:
            if (soln==soln.T).all():
                print('[symmetric found]')
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