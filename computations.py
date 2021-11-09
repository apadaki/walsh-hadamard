import numpy as np

# increment binary list with bits (1,-1)
def increment_binary_list(bin_list):
    idx = len(bin_list) - 1
    while idx >= 0:
        bin_list[idx] = -bin_list[idx]
        if bin_list[idx] == -1:
            break
        idx-=1
    

# compute hamming distance between two binary lists
def hamming_dist(x, y):
    res = 0
    for i in range(len(x)):
        res += (x[i] != y[i])
    return res

# true values of f(H_k), the minimum number of flips to make H_k rank 1
def known(k):
    return [0,1,4,22,96,432,-1,-1,-1][k]

def iterative_comp(H, k):
    # flipset, rowset for finding patterns
    flipset = set()
    rowset = set()
    if k < 1 or k > 7:
        exit(1)
    print('H_{}:'.format(k))
    print(H[k])
    known_min_flips = known(k)
    num_optimal_solns = 0
    min_flips = -1
    dim = 2**k
    
    # initial candidate and empty list of optimal candidates
    candidate = np.ones(dim)
    optimal_candidates = []
    solutions = []
    for i in range(2**(dim)):
        if not i == 0 and i % 1000 == 0:
            print('({}/{})'.format(i,2**dim))
        flips = 0
        # for each row, compute distance between candidate and row (or the inverse candidate if it is closer)
        for j in range(dim):
            dist = hamming_dist(candidate, H[k][j])
            r = min(dist,dim-dist)
            flips += r
        
        # store optimal solutions and count repeats
        if flips == known_min_flips:
            num_optimal_solns+=1
            optimal_candidates.append(candidate.copy())
            
            soln=np.empty((dim,dim),dtype=int)
            for j in range(dim):
                d = hamming_dist(candidate, H[k][j])
                rowset.add(min(d, dim-d))
                inv = False
                if d > dim/2:
                    inv = True
                soln[j] = -candidate.copy() if inv else candidate.copy()
            
            solutions.append(soln)

        flipset.add(flips)
        # check if min_flips must be set to new minimum
        if min_flips == -1 or min_flips > flips:
            min_flips = flips
        
        # move on to next candidate
        increment_binary_list(candidate)

    # print first 32 optimal solutions
    print('\nlist of (<=32) solutions:')
    for soln in solutions[:32]:
        print(np.array(H[k]!=soln, dtype=int), '\n')
    # print set of total flips corresponding all (not necessarily optimal) solutions found by the method of candidate checking. note that values are only checked and stored here if they correspond to an optimal solution GIVEN some candidate.
    print('min_flips: ', min_flips)
    print('num_optimal_solutions:', num_optimal_solns)
    print('flipset:', flipset) 
    print('rowset:', rowset) 
    
def montecarlo_comp(H, k, n_iter=100):
    dim = 2**k
    min_flips = -1

    min_candidate = 0
    num_optimal_solns = 0
    for i in range(n_iter):
        if not i == 0 and i % 200 == 0:
            print('({}/{}) | min_flips: {}'.format(i,n_iter, min_flips))

        # generate random candidate
        candidate = np.random.randint(0,2, size=(dim))*2-1
        flips = 0

        for j in range(dim):
            dist = hamming_dist(candidate, H[k][j])
            r = min(dist,dim-dist)
            flips += r

        if min_flips == -1 or min_flips > flips:
            min_flips = flips
            min_candidate = candidate
            num_optimal_solns = 0
        elif min_flips == flips:
            num_optimal_solns+=1

    print('min_flips: ', min_flips)
    print('proportion_min: {:.10f}'.format(num_optimal_solns/n_iter))
    print('example candidate: ', (min_candidate+1)/2)


if __name__ == '__main__':
    # compute list of H0 through H7
    H = [np.array([1], dtype=int)]
    h1 = np.array([[1,1],[1,-1]])
    for i in range(1,8):
        H.append(np.kron(h1,H[i-1]))

    # accept k >= 1 as imput
    k = int(input('enter value of k (1 to 7): '))

    # iterative_comp(H, k)
    montecarlo_comp(H,k, n_iter = 50000)
