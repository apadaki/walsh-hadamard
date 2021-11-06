#include <iostream>
#include <algorithm>
#include <vector>

#define ull unsigned long long int

using namespace std;


ull two_power(int x) {
    return (ull) 1 << x;
}
void display_vector(vector<ull> vec) {
    for (ull v : vec) {
        printf("%d ", v);
    }
    printf("\n");
}
int num_one_bits(int num) {
    int res = 0;
    while (num > 0) {
        res += (num % 2 == 1);
        num /= 2;
    }
    return res;
}

int min_flips(vector<vector<ull>> H, int k) {
    int dim = two_power(k);
    ull two_dim = two_power(dim);

    int res = -1;
    for (int cand = 0; cand < two_dim; cand++) {
        if (cand % 100000 == 0)
            cout << cand << endl;
        int flips = 0;
        for (int r = 0; r < dim; r++) {
            int d = num_one_bits(cand ^ H[k][r]);
            flips += min(d, dim-d);
        }
        if (res == -1 || res > flips)
            res = flips;
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

    // display_vector(H[5]);

    int k;
    cin >> k;
    cout << min_flips(H, k);

}