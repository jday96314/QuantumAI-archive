#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>

using namespace std;

vector<int> get_qubits(const string& s) {
    vector<int> q;
    if (s.rfind("INFO_TAG", 0) == 0 || s.rfind("REGISTER", 0) == 0 || s.rfind("APPEND_TO_REGISTER", 0) == 0) return q;
    for (size_t i = 0; i < s.size(); i++) {
        if (s[i] == 'q' && (i == 0 || isspace(s[i-1])) && i + 1 < s.size() && isdigit(s[i+1])) {
            size_t j = i + 1;
            int id = 0;
            while (j < s.size() && isdigit(s[j])) {
                id = id * 10 + (s[j] - '0');
                j++;
            }
            q.push_back(id);
            i = j - 1;
        }
    }
    return q;
}

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    ifstream in(argv[1]);
    string line;
    vector<uint64_t> born(20000, UINT64_MAX);
    vector<uint64_t> died(20000, 0);
    uint64_t n = 0;
    int max_qid = 0;

    while (getline(in, line)) {
        for (int q : get_qubits(line)) {
            if (q >= (int)born.size()) {
                born.resize(q + 1000, UINT64_MAX);
                died.resize(q + 1000, 0);
            }
            if (born[q] == UINT64_MAX) born[q] = n;
            died[q] = n;
            if (q > max_qid) max_qid = q;
        }
        n++;
    }

    cout << "Total lines: " << n << endl;
    cout << "Max qubit ID: " << max_qid << endl;

    int pinned = 840;
    int live_at_end = 0;
    for (int q = pinned; q <= max_qid; q++) {
        if (born[q] != UINT64_MAX) {
            if (died[q] > n - 100000) live_at_end++;
        }
    }
    cout << "Ancilla qubits (>= 840) live in last 100k lines: " << live_at_end << " / " << (max_qid - pinned + 1) << endl;

    // Check skyline: calculate number of live wires at regular sample intervals
    int num_samples = 20;
    for (int s = 0; s <= num_samples; s++) {
        uint64_t line_target = (n * s) / num_samples;
        int live = 0;
        for (int q = 0; q <= max_qid; q++) {
            if (born[q] <= line_target && died[q] >= line_target) live++;
        }
        cout << "Progress " << (s * 100 / num_samples) << "% (line " << line_target << "): " << live << " live qubits" << endl;
    }

    return 0;
}
