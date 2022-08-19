#include "experiments.h"

#include <cassert>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <mutex>
#include <queue>
#include <random>
#include <stack>
#include <string>
#include <thread>
#include <unordered_map>

#include "id_queue.h"
#include "point.h"
#include "quad_node.h"
#include "skarf.h"

pair<int, int> Experiments::get_two_random_nodes() {
    size_t n = graph->get_n();
    string seed_str = "very_random_seed";
    seed_seq seed(seed_str.begin(), seed_str.end());
    static mt19937 gen(seed);
    static uniform_int_distribution<> node_distribution(0, n - 1);
    int start_index = node_distribution(gen);
    int destination_index = node_distribution(gen);
    while (start_index == destination_index) {
        destination_index = node_distribution(gen);
    }
    return make_pair(start_index, destination_index);
}

extern bool _compare_with_dijkstra;
void Experiments::run_comparisons(int num_experiments, int partition_size) {
    ofstream file(results_folder + "compare_" + get_graph()->region + "_" + to_string(partition_size) + "_" + to_string(num_experiments) + ".csv");
    assert(file);
    file << "source_id;target_id;source_cell;target_cell;";
    file << "arcflags_path_length;arcflags_visited;arcflags_relaxed;arcflags_querytime;bidirectional_arcflags_path_length;bidirectional_arcflags_visited;bidirectional_arcflags_relaxed;bidirectional_arcflags_querytime;";
    file << "skarf+_path_length;skarf+_visited;skarf+_relaxed;skarf+_querytime;bidirectional_skarf+_path_length;bidirectional_skarf+_visited;bidirectional_skarf+_relaxed;bidirectional_skarf+_querytime\n";

    vector<pair<int, int>> od_pairs(num_experiments);
    vector<string> od_results;
    vector<Edge*> result_skarf;
    vector<vector<Edge*>> dijkstraPaths(num_experiments);
    auto do_query_and_save_result = [&](bool bidirectional, bool arcflags, bool skarf, int query_number) {
        int source_id = od_pairs[query_number].first;
        int destination_id = od_pairs[query_number].second;
        Skarf::get_instance()->reset();
        list<Edge*> temp_result;

        auto startTime = chrono::high_resolution_clock::now();
        if (bidirectional)
            temp_result = Skarf::get_instance()->calculate_route_bidirectional_with(source_id, destination_id, arcflags, skarf);
        else
            temp_result = Skarf::get_instance()->calculate_route_with(source_id, destination_id, arcflags, skarf);
        auto endTime = chrono::high_resolution_clock::now();

        result_skarf.clear();
        result_skarf.reserve(temp_result.size());
        for (Edge* e : temp_result)
            result_skarf.emplace_back(e);

        size_t us = chrono::duration_cast<chrono::microseconds>(endTime - startTime).count();
        if (!dijkstraPaths[query_number].empty() && _compare_with_dijkstra) {
            assert(result_skarf == dijkstraPaths[query_number]);
        }
        if (bidirectional || arcflags || skarf)
            od_results[query_number] += to_string(result_skarf.size()) + ";" + to_string(Skarf::get_instance()->get_visited_edges()) + ";" + to_string(Skarf::get_instance()->get_settled_nodes()) + ";" + to_string(us) + ";";
        else
            dijkstraPaths[query_number] = result_skarf;
    };

    for (auto& od : od_pairs) {
        od = get_two_random_nodes();
        int source_cell_id = graph->get_nodes()[od.first]->cell_idx;
        int target_cell_id = graph->get_nodes()[od.second]->cell_idx;
        od_results.push_back(to_string(od.first) + ";" + to_string(od.second) + ";" + to_string(source_cell_id) + ";" + to_string(target_cell_id) + ";");
    }
    if (_compare_with_dijkstra) {
        for (int i = 0; i < num_experiments; i++) {
            do_query_and_save_result(false, false, false, i);
            assert(dijkstraPaths[i].size());
            if (i % 50 == 0) cout << "Start Dijkstra " << i << " out of " << num_experiments << endl;
        }
    }
    for (int i = 0; i < num_experiments; i++) {
        do_query_and_save_result(false, true, false, i);
        if (i % 50 == 0) cout << "Start Arc-Flags " << i << " out of " << num_experiments << endl;
    }
    for (int i = 0; i < num_experiments; i++) {
        do_query_and_save_result(true, true, false, i);
        if (i % 50 == 0) cout << "Start bidirectional Arc-Flags " << i << " out of " << num_experiments << endl;
    }
    for (int i = 0; i < num_experiments; i++) {
        do_query_and_save_result(false, true, true, i);
        if (i % 50 == 0) cout << "Start SKARF+ " << i << " out of " << num_experiments << endl;
    }
    for (int i = 0; i < num_experiments; i++) {
        do_query_and_save_result(true, true, true, i);
        if (i % 50 == 0) cout << "Start bidirectional SKARF+ " << i << " out of " << num_experiments << endl;
    }
    for (int i = 0; i < num_experiments; i++) od_results[i][od_results[i].size() - 1] = '\n';
    for (auto& l : od_results) file << l;
    file.close();
    return;
}