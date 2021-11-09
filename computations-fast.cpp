#include <bits/stdc++.h>

#define ull unsigned long long int

using namespace std;

struct pair_hash {
    template <class T1, class T2>
    size_t operator () (const pair<T1,T2> &p) const {
        auto h1 = hash<T1>{}(p.first);
        auto h2 = hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;  
    }
};

using Ints = pair<ull, ull>;
using MyUnorderedMap = unordered_map<Ints, int, pair_hash>;


ull two_power(int x) {
    return (ull) 1 << x;
}
void display_vector(vector<ull> vec) {
    for (ull v : vec) {
        printf("%d ", v);
    }
    printf("\n");
}
void display_set(unordered_set<int> myset) {
    for (int s : myset) {
        printf("%d ", s);
    }
    printf("\n");
}
ull num_one_bits(ull num) {
    int res = 0;
    while (num > 0) {
        res += (num % 2 == 1);
        num /= 2;
    }
    return res;
}

// ull num_one_bits_rec(ull num, ull mod_arg, unordered_map<ull,int>& dists_num) {
//     if (dists_num.find(num) != dists_num.end())
//         return dists_num[num];
//     int res = 0;
//     while (num > 0) {
//         res += (num % 2 == 1);
//         num /= 2;

//     }
//     return res;
// }

int hamming_distance(ull a, ull b, ull mod_arg, MyUnorderedMap& dists) {
    pair<ull,ull> p(a,b);
    pair<ull,ull> q(b,a);
    if (mod_arg == 1)
        return a != b;
    if (dists.find(p) == dists.end() && dists.find(q) == dists.end()) {
        int d1 = hamming_distance(a % mod_arg, b % mod_arg, sqrt(mod_arg), dists);
        int d2 = hamming_distance(a / mod_arg, b / mod_arg, sqrt(mod_arg), dists);
        dists[p] = dists[q] = d1+d2;
    }
    if (dists.find(p) == dists.end())
        return dists[q];
    return dists[p];
}

int min_flips_stat(vector<vector<ull>> H, int k, int max_samples) {
    int dim = two_power(k);
    ull two_dim = two_power(dim);
    int res = -1;
    int num_optimal_solutions = 0;
    ull two_32 = ((ull)1 << 32);

    MyUnorderedMap dists;

    for (int i = 0; i < max_samples; i++) {
        if (i != 0 && i % 1000 == 0)
            printf("%d/%d | min_flips: %d | dists_size: %d\n", i, max_samples, res, dists.size());

        int flips = 0;
        ull cand = 0;
        if (k < 6) 
            cand = random() % two_dim;
        else {
            cand = random()*two_power(32) + random();
        }
        for (int r = 0; r < dim; r++) {
            // ull dist = num_one_bits(cand ^ H[k][r]);
            ull dist = hamming_distance(cand, H[k][r], two_32, dists);
            // flips += min(dist, dim-dist);
        }

        if (res == -1 || res > flips) {
            res = flips;
            num_optimal_solutions = 0;
        }
        if (res == flips)
            num_optimal_solutions++;
    }

    cout << num_optimal_solutions << endl;
    printf("proportion_min: {%.8f}\n", (float)num_optimal_solutions/max_samples);

    return res;
}

int min_flips(vector<vector<ull>> H, int k, int known_min, unordered_set<int>& rowset, unordered_set<int>& flipset, int& repeats) {
    int dim = two_power(k);
    ull two_dim = two_power(dim);
    int res = -1;
    repeats = 0;
    for (int cand = 0; cand < two_dim/2; cand++) {
        if (cand != 0 && cand % 1000000 == 0)
            printf("%d/%lld\n", cand, two_dim/2);
        int flips = 0;
        for (int r = 0; r < dim; r++) {
            int dist = num_one_bits(cand ^ H[k][r]);
            flips += min(dist, dim-dist);
        }
        if (res == -1 || res > flips)
            res = flips;
        
        flipset.insert(flips);
        if (flips == known_min) {
            repeats++;
            for (int r = 0; r < dim; r++) {
                int dist = num_one_bits(cand ^ H[k][r]);
                rowset.insert(min(dist, dim-dist));            
            }
        }
    }
    return res;
}
int main() {
    vector<vector<ull>> H(7);
    for (int i = 0; i < 7; i++) {
        int two_i = two_power(i);
        H[i] = vector<ull>(two_i);
        if (i == 0)
            H[i][0] = 0;
        else {
            for (int j = 0; j < two_i/2; j++)
                H[i][j] = H[i-1][j] * (two_power(two_i/2) + 1);
            for (int j = two_i/2; j < two_i; j++)
                H[i][j] = H[i-1][j-two_i/2] * (two_power(two_i/2))-1-H[i-1][j-two_i/2] + two_power(two_i/2);
        }
    }

    // display_vector(H[6]);


    int k;
    cin >> k;
    unordered_set<int> rowset, flipset;
    int repeats = 0;

    // int known_min[] = {0,1,4,22,96,432,-1,-1};
    // printf("min_flips: %d\n", min_flips(H, k, known_min[k], rowset, flipset, repeats));
    // printf("num_optimal_solutions: %d\n", 2*repeats);
    // printf("rowset: ");
    // display_set(rowset);
    // printf("flipset: ");
    // display_set(flipset);

    printf("min_flips: %d\n", min_flips_stat(H, k, 1000000));

    
}
