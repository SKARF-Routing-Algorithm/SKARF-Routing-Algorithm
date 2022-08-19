#pragma once
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <mutex>
#include <thread_pool.hpp>

#include "graph.h"
#include "id_queue.h"

using namespace std;
using BitLabel = boost::dynamic_bitset<>;

class Skarf {
   private:
    static Skarf* _instance;
    static mutex _mutex;
    Skarf() {}
    ~Skarf() {}

    unsigned visited_edges;
    unsigned settled_nodes;

    vector<BitLabel> arcflags_key_to_label;
    vector<BitLabel> skarf_key_to_label;
    vector<BitLabel> preprocessing_key_to_label;
    vector<size_t> arcflags_edge_to_key;
    vector<size_t> skarf_edge_to_key;
    vector<size_t> preprocessing_edge_to_key;

    static const int num_steps = 20000;
    static long long step_size;
    void test();
    unsigned num_threads = 1;
    int partition_size;
    bool precomputed = false;

    bool is_set_skarf(int edge, int startCell, int destinationCell, int flagType);
    bool is_set_skarf_plus(int edge, int startCell, int destinationCell);
    bool is_set_arcflags(int edge, int destinationCell);
    bool is_set_skarf_bidirectional(int edge, int startCell, int destinationCell, int flagType, bool forwardSearch);
    bool is_set_skarf_plus_forward(int edge, int startCell, int destinationCell);
    bool is_set_skarf_plus_backward(int edge, int startCell, int destinationCell);
    bool is_set_arcflags_forward(int edge, int startCell, int destinationCell);
    bool is_set_arcflags_backward(int edge, int startCell, int destinationCell);
    static void precompute_cell(int cell_idx, int partition_size);
    static void precompute_direction(int cell_idx, bool transposed, int partition_size, unordered_map<Node*, long long>& boundary_dist);
    static void save_flags(int cell_idx, string folder, BitLabel& flagged, int partition_size);
    static void mark_edges(Node* source, BitLabel& skeleton_union, BitLabel& arc_flags_union, bool transposed, unordered_map<Node*, long long>& boundary_dist);
    static void add_cell_overlap(int cell_idx, bool transposed, int partition_size, Node* center);
    void export_overlap_results();
    void set_step_size();

    static Graph* graph;
    static SkarfQueue::MinIDQueue forward_queue;
    static SkarfQueue::MinIDQueue backward_queue;

   public:
    Skarf(Skarf& to_clone) = delete;
    void operator=(const Skarf&) = delete;
    static Skarf* get_instance();
    static string results_folder;
    void set_partition_size(int size) { partition_size = size; };
    void merge_flags(string folder);
    void import_flags(int flagType);
    void export_flags(string folder);
    void precompute(int start, int end);

    static void set_graph(Graph* g);
    static Graph* get_graph() { return graph; };

    unsigned get_settled_nodes() { return settled_nodes; };
    unsigned get_visited_edges() { return visited_edges; };
    void settle_node() { ++settled_nodes; };
    void visit_edge() { ++visited_edges; };

    list<Edge*> calculate_route_with(int source_id, int target_id, bool arcflags, bool skarf);
    list<Edge*> calculate_route_bidirectional_with(int source_id, int target_id, bool arcflags, bool skarf);
    void compress(unordered_map<Edge*, BitLabel>& labels);
    void import_all_flags(bool import_arc_flags = true, bool import_skarf_flags = true);
    void backup_flags(string flag_type);
    void set_threads(unsigned num) { num_threads = num; };
    void reset();
};