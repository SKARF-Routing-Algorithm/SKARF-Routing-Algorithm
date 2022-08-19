#pragma once
#include <vector>
using namespace std;

struct SCC_Graph {
    vector<vector<int>> adj;
    vector<int> scc_index;
    void compute_scc_graph(vector<vector<int>>& graph);
};
