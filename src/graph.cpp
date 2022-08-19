#include "graph.h"

#include <id_queue.h>

#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "point.h"
#include "quad_node.h"
#include "scc.h"
#include "util.h"
#include "metis.h"
#include "thread_pool.hpp"

using namespace std;

Graph::Graph() {
    top_left = new Point(180, -90);
    bot_right = new Point(-180, 90);
}

void Graph::update_min_max_lat(double lat) {
    if (lat > top_left->get_lat())
        top_left->set_lat(lat);
    if (lat < bot_right->get_lat())
        bot_right->set_lat(lat);
}

void Graph::update_min_max_lon(double lon) {
    if (lon < top_left->get_lon())
        top_left->set_lon(lon);
    if (lon > bot_right->get_lon())
        bot_right->set_lon(lon);
}

void Graph::load_network(string region_folder, bool preprocess_raw, unsigned num_threads) {
    this->region_folder = region_folder;
    cout << "loading street network..." << endl;
    ifstream file(region_folder + "/edges_preprocessed.csv");
    if (!file || preprocess_raw) {
        cout << "preprocessing raw network data..." << endl;
        preprocess_raw = true;
        file.open(region_folder + "/edges_raw.csv");
    }
    string line;
    double current_id = 0;
    getline(file, line);  // skip first line

    map<Point, int> point_index;
    auto get_node = [&](string lon, string lat) -> Node* {
        Point p(stod(lon), stod(lat));
        if (point_index.find(p) == point_index.end()) {
            point_index[p] = nodes.size();
            Node* node = new Node(nodes.size(), p.get_lat(), p.get_lon());
            node->s_lon = lon;
            node->s_lat = lat;
            node->id = nodes.size();
            nodes.push_back(node);
        }
        return nodes[point_index[p]];
    };

    map<pair<Node*, Node*>, Edge*> index_to_edge;
    auto create_edge = [&](Node* a, Node* b, long long time_in_ns, double km, double kmh) {
        if (a == b) return;
        if (index_to_edge.find({a, b}) == index_to_edge.end()) {
            Edge* edge = new Edge(a, b);
            edge->cost = time_in_ns;
            edge->km = km;
            edge->kmh = kmh;
            a->adj.push_back(edge);
            b->back_adj.push_back(edge);
            ++(a->degree);
            ++(b->degree);
            edge->id = edges.size();
            edges.push_back(edge);
            index_to_edge[{a, b}] = edge;
        } else {
            auto [_, edge] = *(index_to_edge.find({a, b}));
            edge->cost = min(edge->cost, time_in_ns);
        }
    };

    while (getline(file, line)) {
        vector<string> csv = split(line, ",", false);
        if (preprocess_raw) {
            double cost = stod(csv[static_cast<int>(RAW_EDGE::cost)]);
            double reverse_cost = stod(csv[static_cast<int>(RAW_EDGE::reverse_cost)]);
            double km = stod(csv[static_cast<int>(RAW_EDGE::km)]);
            double kmh = stod(csv[static_cast<int>(RAW_EDGE::kmh)]);
            string x1 = csv[static_cast<int>(RAW_EDGE::x1)];
            string y1 = csv[static_cast<int>(RAW_EDGE::y1)];
            string x2 = csv[static_cast<int>(RAW_EDGE::x2)];
            string y2 = csv[static_cast<int>(RAW_EDGE::y2)];
            Node* from = get_node(x1, y1);
            Node* to = get_node(x2, y2);
            long long time_in_ns = (long long)(cost * 3600 * 1e9);
            create_edge(from, to, time_in_ns, km, kmh);
            if (reverse_cost <= 1e5) {
                long long reverse_time_in_ns = (long long)(cost * 3600 * 1e9);
                create_edge(to, from, reverse_time_in_ns, km, kmh);
            }
        } else {
            string x1 = csv[static_cast<int>(PREPROCESSED_EDGE::x1)];
            string y1 = csv[static_cast<int>(PREPROCESSED_EDGE::y1)];
            string x2 = csv[static_cast<int>(PREPROCESSED_EDGE::x2)];
            string y2 = csv[static_cast<int>(PREPROCESSED_EDGE::y2)];
            double km = stod(csv[static_cast<int>(PREPROCESSED_EDGE::km)]);
            double kmh = stod(csv[static_cast<int>(PREPROCESSED_EDGE::kmh)]);
            double cost = stod(csv[static_cast<int>(PREPROCESSED_EDGE::cost)]);
            Node* from = get_node(x1, y1);
            Node* to = get_node(x2, y2);
            create_edge(from, to, cost, km, kmh);
        }
    }
    file.close();
    if (preprocess_raw) {
        cout << "before preprocessing: " << get_n() << " nodes, " << get_m() << " edges" << endl;
        filter_main_component();
        filter_shortest_edges(num_threads);
        perturb_edge_costs();
        write_preprocessed(region_folder);
        verify_uniqueness_of_shortest_paths_between_boundary_nodes(num_threads);
    }
    verify_edge_consistency();
    cout << "preprocessed graph has " << get_n() << " nodes and " << get_m() << " edges" << endl;
    build_search_tree();
}

void Graph::write_preprocessed(string region_folder) {
    set<Point> points;
    for (Node* node : nodes) {
        Point p(node->lon, node->lat);
        assert(points.find(p) == points.end());
        points.insert(p);
    }

    cout << "writing preprocessed network data..." << endl;
    ofstream out(region_folder + "/edges_preprocessed.csv");
    out << "x1,y1,x2,y2,km,kmh,cost_in_ns" << endl;
    for (Edge* e : edges) {
        out << e->from->s_lon << "," << e->from->s_lat << ",";
        out << e->to->s_lon << "," << e->to->s_lat << ",";
        out << e->km << "," << e->kmh << "," << e->cost << endl;
    }
    out.close();
}

void Graph::partition(int nparts) {
    // METIS expects an undirected graph.
    // So we input the underlying undirected graph.
    set<pair<Node*, Node*>> arc_present;
    for (Edge* edge : edges)
        arc_present.insert({edge->from, edge->to});
    vector<vector<int>> undirected_adj(nodes.size());
    for (Edge* edge : edges) {
        undirected_adj[edge->from->id].push_back(edge->to->id);
        auto reverse_arc = make_pair(edge->to, edge->from);
        if (arc_present.find(reverse_arc) == arc_present.end())
            undirected_adj[edge->to->id].push_back(edge->from->id);
    }

    vector<int> xadj(get_n() + 1), adjncy;
    xadj[0] = 0;

    for (int i = 0; i < int(nodes.size()); ++i) {
        for (int node_id : undirected_adj[i])
            adjncy.push_back(node_id);
        xadj[i + 1] = adjncy.size();
    }

    vector<int> partitions;
    partitions.resize(get_n());
    cout << "Graph formatted" << endl;

    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;

    int wgtflag = 0, numflag = 1, edgecut = 0;
    int ncon = 1;
    auto r = &(partitions[0]);
    int tempN = get_n();
    METIS_PartGraphRecursive(&tempN, &ncon, &(xadj[0]), &(adjncy[0]), NULL, NULL, NULL, &nparts, NULL, NULL, options, &edgecut, &(partitions[0]));
    cout << "finished" << endl;

    parts.clear();
    parts.resize(nparts);
    boundary_nodes.clear();
    boundary_nodes.resize(parts.size());

    for (int i = 0; i < int(partitions.size()); i++) {
        int idx = partitions[i];
        assert(idx < nparts);
        parts[idx].push_back(nodes[i]);
        nodes[i]->cell_idx = idx;
    }
    gather_boundary_nodes();
}

void Graph::gather_boundary_nodes() {
    cout << "compute boundary nodes" << endl;
    unordered_map<long long, int> partMap;

    for (int i = 0; i < parts.size(); i++)
        for (auto node : parts[i])
            partMap[node->id] = i;

    int maxCount = 0;
    for (int i = 0; i < parts.size(); i++) {
        int count = 0;
        for (auto node : parts[i]) {
            for (auto edge : node->adj) {
                if (partMap[edge->to->id] != i) {
                    node->boundaryNode = true;
                    boundary_nodes.at(i).push_back(node);
                    count++;
                    break;
                }
            }
            if (!node->boundaryNode)
                for (auto edge : node->back_adj) {
                    if (partMap[edge->from->id] != i) {
                        node->boundaryNode = true;
                        boundary_nodes.at(i).push_back(node);
                        count++;
                        break;
                    }
                }
        }
        if (count > maxCount)
            maxCount = count;
    }
    max_boundary_nodes = maxCount;
    cout << "finished computing boundary nodes" << endl;
}

void Graph::build_search_tree() {
    if (get_n() == 0 && get_m() == 0) return;
    root_node = new QuadNode(top_left, bot_right);
    cout << "building searchtree!" << endl;
    for (Node* node : nodes) {
        if (node->lat && node->lon)
            root_node->insert(node);
    }
    cout << "searchtree successfully build!" << endl;
}

void Graph::get_max_degree() {
    max_degree = 0;
    for (Node* node : nodes)
        max_degree = max(max_degree, node->degree);
}

string Graph::get_current_time_point_as_string() {
    time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now());
    string ts = ctime(&t);
    ts.resize(ts.size() - 1);
    return ts;
}

void Graph::export_region_to_csv(int index) {
    string time = get_current_time_point_as_string();
    ofstream file("region_ " + to_string(index) + "_" + time + ".csv");
    file << "lat;lon;cell_idx;deg\n";
    for (Node* node : parts[index]) {
        file << node->lat << ";" << node->lon << ";";
        file << index << ";" << node->degree << "\n";
    }
    file.close();
}

void Graph::export_partition_to_csv() {
    ofstream file("region.csv");
    file << "lat;lon;cell_idx;deg\n";
    for (int index = 0; index < (int)parts.size(); ++index)
        for (Node* node : parts[index]) {
            file << node->lat << ";" << node->lon << ";";
            file << index << ";" << node->degree << "\n";
        }
    file.close();
}

void Graph::export_edge_hit_count_to_csv(unordered_map<Edge*, int>& hit_count, int threshold = 1) {
    string time = get_current_time_point_as_string();
    ofstream file("edges_" + time + ".csv");
    file << "from_lat;from_lon;to_lat;to_lon;value\n";
    for (auto id : node_ids)
        for (auto edge : nodes[id]->adj)
            if (hit_count[edge] >= threshold) {
                file << edge->from->lat << ";" << edge->from->lon << ";";
                file << edge->to->lat << ";" << edge->to->lon << ";";
                file << hit_count[edge] << "\n";
            }
    file.close();
}

void Graph::export_edge_ratio_to_csv(unordered_map<Edge*, double>& ratio, double threshold) {
    string time = get_current_time_point_as_string();
    ofstream file("edges_" + time + ".csv");
    file << "from_lat;from_lon;to_lat;to_lon;value\n";
    for (auto edge : edges)
        if (ratio[edge] >= threshold) {
            double r = ratio[edge] * 1000;
            file << edge->from->lat << ";" << edge->from->lon << ";";
            file << edge->to->lat << ";" << edge->to->lon << ";";
            if (r >= 999999)
                file << 999999 << "\n";
            else if (r <= 0.00001)
                file << 0.00001 << "\n";
            else
                file << r << "\n";
        }
    file.close();
}

struct Node* Graph::find(Point p, QuadNode* root) {
    Point* tl = root->get_top_left();
    Point* br = root->get_bot_right();
    double y1 = tl->get_lat();
    double y2 = br->get_lat();
    double x1 = tl->get_lon();
    double x2 = br->get_lon();

    double dist = root->euclidean_distance(tl->get_lat(), tl->get_lon(), br->get_lat(), br->get_lon());
    Result best = Result(NULL, dist);
    best = root->nearest(&p, best);
    return best.node;
}

int majority_element(vector<int>& a) {
    int element = -1;
    int counter = 0;
    for (int e : a)
        if (e == element)
            ++counter;
        else if (counter == 0) {
            element = e;
            counter = 1;
        } else
            --counter;
    return element;
}

void Graph::filter_main_component() {
    cout << "filter biggest strongly connected component..." << endl;
    vector<vector<int>> pure_index_graph(nodes.size());
    for (Node* v : nodes)
        for (Edge* e : v->adj) {
            pure_index_graph[v->id].push_back(e->to->id);
        }
    SCC_Graph scc;
    scc.compute_scc_graph(pure_index_graph);
    int main_scc = majority_element(scc.scc_index);

    vector<Node*> main_nodes;
    for (Node* v : nodes)
        if (scc.scc_index[v->id] == main_scc) {
            main_nodes.push_back(v);
        }
    edges.clear();
    for (Node* v : main_nodes) {
        vector<Edge*> adj;
        for (Edge* e : v->adj)
            if (scc.scc_index[e->to->id] == main_scc) {
                adj.push_back(e);
                e->id = edges.size();
                edges.push_back(e);
            }
        v->adj = adj;
        assert(adj.size());

        vector<Edge*> back_adj;
        for (Edge* e : v->back_adj)
            if (scc.scc_index[e->from->id] == main_scc) {
                back_adj.push_back(e);
            }
        v->back_adj = back_adj;
        assert(back_adj.size());
    }
    nodes = main_nodes;

    for (int i = 0; i < int(nodes.size()); ++i) {
        nodes[i]->id = i;
    }

    // update edges vector
    vector<Edge*> main_edges;
    for (Node* node : nodes)
        for (Edge* e : node->adj)
            main_edges.push_back(e);
    for (int i = 0; i < int(main_edges.size()); ++i)
        main_edges[i]->id = i;

    node_ids.resize(nodes.size());
    iota(node_ids.begin(), node_ids.end(), 0);
}

void Graph::filter_shortest_edges(unsigned num_threads) {
    cout << "filter shortest edges" << endl;
    thread_pool pool(num_threads);
    // removes all arcs (a,b) which are not the shortest path from a to b

    auto loop = [&](const int begin, const int end) {
        for (int i = begin; i < end; ++i) {
            Node* source = nodes[i];
            int unsettled_neighbors = source->adj.size();
            set<Node*> neighbors;
            for (Edge* e : source->adj) neighbors.insert(e->to);

            int source_id = source->id;
            unordered_map<Node*, long long> dist;
            unordered_map<Node*, Node*> pred;
            unordered_map<Node*, bool> visited;

            SkarfQueue::MinIDQueue queue(nodes.size());
            dist[source] = 0;
            queue.push({source_id, 0});

            while (!queue.empty() && unsettled_neighbors) {
                auto [vi, v_dist] = queue.pop();
                Node* v = nodes[vi];
                if (neighbors.find(v) != neighbors.end()) --unsettled_neighbors;
                for (Edge* e : v->adj) {
                    Node* w = e->to;
                    if (visited[w]) continue;
                    if (dist[w] == 0 || v_dist + e->cost < dist[w]) {
                        long long w_dist = v_dist + e->cost;
                        dist[w] = w_dist;
                        pred[w] = v;
                        if (queue.contains_id(w->id))
                            queue.decrease_key({w->id, w_dist});
                        else
                            queue.push({w->id, w_dist});
                    }
                }
                visited[v] = true;
            }

            vector<Edge*> main_adj;
            for (Edge* e : source->adj)
                if (pred[e->to] == source) {
                    main_adj.push_back(e);
                }
            source->adj = main_adj;
        }
    };
    pool.parallelize_loop(0, nodes.size(), loop);
    pool.wait_for_tasks();
    
    vector<Edge*> main_edges;
    // update back adjacency lists
    for (Node* node : nodes) {
        node->back_adj.clear();
        for (Edge* e : node->adj) {
            main_edges.push_back(e);
        }
    }
    for (Node* from : nodes)
        for (Edge* e : from->adj)
            e->to->back_adj.push_back(e);

    // update edges vector
    edges = main_edges;
    for (int i = 0; i < int(edges.size()); ++i)
        edges[i]->id = i;
}

void Graph::perturb_edge_costs() {
    cout << "perturbing edge costs..." << endl;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> epsilon(0, 1e9);
    for (auto edge : edges) edge->cost += epsilon(gen);
}

void Graph::verify_uniqueness_of_shortest_paths_between_boundary_nodes(unsigned num_threads) {
    cout << "verifying uniqueness of shortest paths between boundary nodes..." << endl;
    thread_pool pool(num_threads);

    auto verify_cell = [&](int cellId) {
        auto& cell_boundary_nodes = boundary_nodes[cellId];
        for (Node* source : cell_boundary_nodes) {
            // compute an arbitrary shortest path tree from source
            vector<long long> dist(nodes.size(), numeric_limits<long long>::max());
            vector<int> pred(nodes.size());
            vector<Edge*> predEdge(nodes.size());
            vector<bool> visited(nodes.size());
            typedef pair<long long, int> queue_element;
            auto cmp = [](queue_element& a, queue_element& b) {
                return a.first > b.first;
            };
            priority_queue<queue_element, vector<queue_element>, decltype(cmp)> q(cmp);
            int src = source->id;
            dist[src] = 0;
            q.emplace(0, src);
            while (!q.empty()) {
                auto [v_dist, vi] = q.top();
                q.pop();
                if (visited[vi]) continue;
                assert(v_dist == dist[vi]);
                if (vi != src) assert(dist[vi] == dist[pred[vi]] + predEdge[vi]->cost);
                for (Edge* e : nodes[vi]->adj) {
                    int wi = e->to->id;
                    if (visited[wi]) continue;
                    if (v_dist + e->cost < dist[wi]) {
                        dist[wi] = v_dist + e->cost;
                        pred[wi] = vi;
                        predEdge[wi] = e;
                        q.emplace(dist[wi], wi);
                    }
                }
                visited[vi] = true;
            }
            // check that the shortest path tree is unique
            for (int vi = 0; vi < get_n(); ++vi) {
                assert(visited[vi]);
                if (vi == src) continue;
                // check if v has a unique predecessor
                int count_optimal_predecessors = 0;
                for (Edge* e : nodes[vi]->back_adj) {
                    int wi = e->from->id;
                    if (dist[wi] + e->cost == dist[vi])
                        ++count_optimal_predecessors;
                }
                assert(count_optimal_predecessors >= 1);
                if (count_optimal_predecessors > 1) {
                    cout << "cells checked: " << cellId + 1 << endl;
                    assert(false);
                }
            }
        }
    };

    auto start = chrono::high_resolution_clock::now();
    for (int cellId = 0; cellId < int(boundary_nodes.size()); ++cellId) {
        pool.push_task(verify_cell, cellId);
    }
    synced_stream sync_out;
    while (pool.get_tasks_total() > 0) {
        sync_out.println(pool.get_tasks_total(),
                         " tasks total, ",
                         pool.get_tasks_running(),
                         " tasks running, ",
                         pool.get_tasks_queued(),
                         " tasks queued.");
        this_thread::sleep_for(std::chrono::milliseconds(5000));
    }
    pool.wait_for_tasks();
    auto end = chrono::high_resolution_clock::now();
    auto minutes = chrono::duration_cast<chrono::minutes>(end - start);
    auto seconds = chrono::duration_cast<chrono::seconds>(end - start);
    cout << "verification of shortest path uniqueness between boundary nodes took ";
    cout << minutes.count() << " minutes and " << seconds.count() << " seconds" << endl;
}

void Graph::verify_edge_consistency() {
    cout << "Verify edge consistency" << endl;
    for (int i = 0; i < get_n(); i++) {
        Node* v = nodes[i];
        for (auto& e_to : nodes[i]->adj) {
            Node* w = e_to->other(v);
            bool back_edge_exists = false;
            for (auto& e_back : w->back_adj) {
                if (e_back->id == e_to->id) back_edge_exists = true;
            }
            assert(back_edge_exists);
        }
        for (auto& e_back : nodes[i]->back_adj) {
            Node* w = e_back->other(v);
            bool to_edge_exists = false;
            for (auto& e_to : w->adj) {
                if (e_to->id == e_back->id) to_edge_exists = true;
            }
            assert(to_edge_exists);
        }
    }
    cout << "All edges are consistent" << endl;
}