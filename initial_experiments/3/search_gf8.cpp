// Exhaustive rank-five bilinear decomposition of GF(8) multiplication.
#include <array>
#include <iostream>

int gf_mul(int a, int b) {
    int p = 0;
    for (int i = 0; i < 3; ++i) if ((b >> i) & 1) p ^= a << i;
    for (int i = 4; i >= 3; --i) if ((p >> i) & 1) p ^= 0b1011 << (i - 3);
    return p;
}

int main() {
    std::array<int, 49> terms{};
    std::array<int, 3> target{};
    for (int k = 0; k < 3; ++k)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                target[k] |= ((gf_mul(1 << i, 1 << j) >> k) & 1) << (3 * i + j);
    for (int a = 1; a < 8; ++a)
        for (int b = 1; b < 8; ++b)
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    terms[(a - 1) * 7 + b - 1] |= (((a >> i) & (b >> j) & 1) << (3 * i + j));
    for (int a = 0; a < 45; ++a)
    for (int b = a + 1; b < 46; ++b)
    for (int c = b + 1; c < 47; ++c)
    for (int d = c + 1; d < 48; ++d)
    for (int e = d + 1; e < 49; ++e) {
        int ids[5] = {a,b,c,d,e};
        int masks[3] = {};
        bool ok = true;
        for (int k = 0; k < 3 && ok; ++k) {
            bool found = false;
            for (int mask = 0; mask < 32 && !found; ++mask) {
                int v = 0;
                for (int q = 0; q < 5; ++q) if ((mask >> q) & 1) v ^= terms[ids[q]];
                if (v == target[k]) found = true, masks[k] = mask;
            }
            ok = found;
        }
        if (!ok) continue;
        for (int q = 0; q < 5; ++q) {
            int id = ids[q];
            int out = ((masks[0] >> q) & 1) | (((masks[1] >> q) & 1) << 1) |
                      (((masks[2] >> q) & 1) << 2);
            std::cout << '{' << id / 7 + 1 << ',' << id % 7 + 1 << ',' << out << "},\n";
        }
        return 0;
    }
    std::cout << "none\n";
}
