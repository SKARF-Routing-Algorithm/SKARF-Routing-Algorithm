#pragma once

#include <queue>
#include <vector>

#include "graph.h"

class Experiments {
   private:
    Graph* graph;

   public:
    string results_folder;
    Experiments(Graph* _g) : graph{_g} {};
    Graph* get_graph() { return graph; };
    pair<int, int> get_two_random_nodes();
    void run_comparisons(int num_experiments, int partition_size);
};
