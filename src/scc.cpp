#include "scc.h"
#include<functional>
#include<numeric>
#include <cassert>
#include <sys/resource.h>
#include <iostream>

void ensure_sufficient_stack_size() {
    rlim_t stack_size = 1'000'000'000;
    rlimit limit;
    int answer = getrlimit(RLIMIT_STACK, &limit);
    if (!answer) {
        if (limit.rlim_cur < stack_size) {
            limit.rlim_cur = stack_size;
            answer = setrlimit(RLIMIT_STACK, &limit);
            if (answer != 0) {
                std::cerr << "increasing the stack size went wrong, error code: " << answer << std::endl;
            }
        }
    } else {
        std::cerr << "getting the current current stack size failed" << std::endl;
    }
}

void SCC_Graph::compute_scc_graph(vector<vector<int>>& graph) {
	ensure_sufficient_stack_size();

	int n = int(graph.size());
	// get finishing times
	vector<bool> visited(n, false);
	vector<int> f(n);
	int time = 0;
	function<void(int)> dfs = [&](int v) {
		++time;
		assert(v < n);
		visited[v] = true;
		for (int w : graph[v]) {
			assert(w < n);
			if (!visited[w]) dfs(w);
		}
		++time;
		f[v] = time;
	};
	for (int v = 0; v < n; ++v)
		if (!visited[v]) dfs(v);
	// collect components
	vector<int> order(n);
	iota(order.begin(), order.end(), 0);
	sort(order.begin(), order.end(), [&](int i, int j) {
		return f[i] > f[j];
	});
	visited = vector<bool>(n, false);
	int scc = 0;
	scc_index.resize(n);
	vector<vector<int>> transposed_graph(n);
	for (int v = 0; v < n; ++v) for (int w : graph[v])
		transposed_graph[w].push_back(v);
	function<void(int)> collect_components = [&](int v) {
		visited[v] = true;
		scc_index[v] = scc;
		for (int w : transposed_graph[v])
			if (!visited[w]) collect_components(w);
	};
	for (int v : order) {
		if (!visited[v]) {
			collect_components(v);
			++scc;
		}
	}
	// generate arcs between components
	adj.resize(scc);
	for (int v = 0; v < n; ++v)
		for (int w : graph[v])
			if (scc_index[v] != scc_index[w])
				adj[scc_index[v]].push_back(scc_index[w]);
	// clean up redundant arcs
	for (int comp = 0; comp < scc; ++comp) {
		sort(adj[comp].begin(), adj[comp].end());
		auto last = unique(adj[comp].begin(), adj[comp].end());
		adj[comp].erase(last, adj[comp].end());
	}
}