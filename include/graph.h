#pragma once
#include <string>
#include <unordered_map>
#include "quad_node.h"
#include <bitset>
#include <set>
using namespace std;

enum class RAW_EDGE { id, osm_source_id, osm_target_id, km, kmh, cost, reverse_cost, x1, y1, x2, y2 };
enum class PREPROCESSED_EDGE { x1, y1, x2, y2, km, kmh, cost };

using namespace std;

class Graph {
private:
    int max_degree;
    QuadNode* root_node;
    Point* top_left;
    Point* bot_right;
    vector<Node*> nodes;
    vector<Edge*> edges;
    vector<int> node_ids;
    vector<vector<Node*>> parts;
    vector<vector<Node*>> boundary_nodes;
    int max_boundary_nodes;
    void update_min_max_lat(double lat);
    void update_min_max_lon(double lon);
    void gather_boundary_nodes();
    void get_max_degree();
    Node* unvisited_bodes(unordered_map<long long, bool> nodes);
    void filter_shortest_edges(unsigned num_threads);
public:
    string region_folder;
    string region;
    Graph();
    struct Edge* create_edge(struct Node* from, struct Node* to, double km, double kmh, double cost, long long id);
    void read_edges(string path);
    void load_network(string edgePath, bool preprocess_raw = false, unsigned num_threads = 1);
    struct Node* find(Point p, QuadNode* root);
    void build_search_tree();
    void partition(int nparts);
    void filter_main_component();
    string get_current_time_point_as_string();
    void perturb_edge_costs();
    void verify_uniqueness_of_shortest_paths_between_boundary_nodes(unsigned num_threads);
    void verify_edge_consistency();

    void write_preprocessed(string region_folder);
    void export_region_to_csv(int index);
    void export_edge_hit_count_to_csv(unordered_map<Edge*, int>& hit_count, int threshold);
    void export_partition_to_csv();
    void export_arc_flags_edge_marker();
    void export_edges_to_gr(unordered_map<Edge*, bool>& hit_count, int cell_idx);
    void export_edge_ratio_to_csv(unordered_map<Edge*, double>& ratio, double threshold = 0);

    int get_n() { return int(nodes.size()); }
    size_t get_m() { return edges.size(); }
    int get_max_boundary_nodes() { return max_boundary_nodes; }
    vector<int>& get_node_ids() { return node_ids; }
    vector<vector<Node*>>& get_parts() { return parts; }
    vector<Node*>& get_nodes() { return nodes; }
    vector<Edge*>& get_edges() { return edges; }
    QuadNode* get_street_node_root_node() { return root_node; }
    vector<vector<Node*>>& get_boundary_node_lists() { return boundary_nodes; }
};