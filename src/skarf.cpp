#include "skarf.h"

#include <pthread.h>
#include <sys/types.h>
#include <unistd.h>

#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <mutex>
#include <stack>
#include <stdexcept>
#include <string>
#include <thread>

#include "graph.h"
#include "id_queue.h"
#include "segment_tree.h"
#include "util.h"
using namespace std;

Skarf* Skarf::_instance = nullptr;
mutex Skarf::_mutex;
long long Skarf::step_size;
string Skarf::results_folder;

Skarf* Skarf::get_instance() {
    lock_guard<mutex> lock(_mutex);
    if (_instance == nullptr) _instance = new Skarf();
    return _instance;
}

void Skarf::mark_edges(Node* source, BitLabel& skeleton_union, BitLabel& arc_flags_union, bool transposed, unordered_map<Node*, long long>& boundary_dist) {
    int N = graph->get_n();
    auto& nodes = graph->get_nodes();

    vector<long long> dist(N, numeric_limits<long long>::max());
    vector<long long> reach(N);
    vector<int> pred(N);
    vector<Edge*> pred_edge(N);
    vector<bool> visited(N);
    vector<vector<int>> children(N);

    SkarfQueue::MinIDQueue queue(N);
    int src = source->id;
    dist[src] = 0;
    queue.push({src, 0});

    while (!queue.empty()) {
        auto [vi, v_dist] = queue.pop();
        assert(v_dist == dist[vi]);
        if (vi != src) {
            assert(dist[vi] == dist[pred[vi]] + pred_edge[vi]->cost);
            children[pred[vi]].push_back(vi);
            arc_flags_union[pred_edge[vi]->id] = 1;
        }
        auto& neighbors = transposed ? nodes[vi]->back_adj : nodes[vi]->adj;
        for (Edge* e : neighbors) {
            int wi = e->other(vi);
            if (visited[wi]) continue;
            if (v_dist + e->cost < dist[wi]) {
                dist[wi] = v_dist + e->cost;
                pred[wi] = vi;
                pred_edge[wi] = e;
                if (queue.contains_id(wi))
                    queue.decrease_key({wi, dist[wi]});
                else
                    queue.push({wi, dist[wi]});
            }
        }
        visited[vi] = true;
    }

    for (auto [node, _] : boundary_dist)
        boundary_dist[node] = max(boundary_dist[node], dist[node->id]);

    vector<bool> visited_reach_computation(N);
    stack<int> s;
    s.emplace(src);
    while (!s.empty()) {
        int vi = s.top();
        if (!visited_reach_computation[vi]) {
            assert(reach[vi] == 0);
            for (int wi : children[vi]) {
                s.emplace(wi);
            }
            visited_reach_computation[vi] = true;
        } else {
            s.pop();
            // set the flags
            for (int wi : children[vi]) {
                Edge* e = pred_edge[wi];
                assert(reach[wi] + e->cost > 0);
                assert(dist[wi] == dist[vi] + e->cost);
                reach[vi] = max(reach[vi], e->cost + reach[wi]);
                // see also our skeleton definition for non-geometric-realisations
                if ((reach[wi] + e->cost) > dist[vi]) {
                    skeleton_union[e->id] = 1;
                }
            }
        }
    }
}

void Skarf::save_flags(int cell_idx, string folder, BitLabel& flagged, int partition_size) {
    string dir = graph->region_folder + "/intermediate/" + folder;
    validate_directory(dir);
    ofstream out(dir + "/flags_partitionSize" + to_string(partition_size) + "_cell" + to_string(cell_idx));
    for (int i = 0; i < int(graph->get_m()); ++i)
        if (flagged[i]) out << i << endl;
    out.close();
};

void Skarf::precompute_direction(int cell_idx, bool transposed, int partition_size, unordered_map<Node*, long long>& boundary_dist) {
    BitLabel skeleton_union(graph->get_m(), 0);
    BitLabel arc_flags_union(graph->get_m(), 0);

    auto& cell = graph->get_parts()[cell_idx];
    for (auto node : cell) {
        if (node->boundaryNode) {
            mark_edges(node, skeleton_union, arc_flags_union, transposed, boundary_dist);
        }
        for (auto edge : node->adj) {
            if (edge->to->cell_idx == cell_idx) {
                skeleton_union[edge->id] = 1;
                arc_flags_union[edge->id] = 1;
            }
        }
    }
    int index = cell_idx;
    if (!transposed) index += partition_size;
    save_flags(index, "arcflags", arc_flags_union, partition_size);
    save_flags(index, "skarf", skeleton_union, partition_size);
}

int maximum_overlap(vector<pair<long long, long long>>& intervals) {
    // treating intervals as open on the left side and closed on the right side
    vector<pair<long long, int>> events;
    for (auto [l, r] : intervals) {
        events.emplace_back(l, 1);  // open event
        events.emplace_back(r, 0);  // close event
    }
    sort(events.begin(), events.end());
    int maximum = 0;
    int counter = 0;
    for (auto [_, type] : events)
        if (type)
            ++counter;
        else {
            maximum = max(maximum, counter);
            --counter;
        }
    return maximum;
}

void Skarf::add_cell_overlap(int cell_idx, bool transposed, int partition_size, Node* center) {
    string file = graph->region_folder + "/intermediate/skarf/flags_partitionSize" + to_string(partition_size) + "_cell" + to_string(cell_idx);
    ifstream in(file);
    vector<int> skeleton_union;
    string line;
    while (getline(in, line)) {
        skeleton_union.push_back(atoi(line.c_str()));
    }
    in.close();

    // compute shortest path tree from center node
    int N = graph->get_n();
    auto& nodes = graph->get_nodes();

    vector<long long> dist(N, numeric_limits<long long>::max());
    vector<long long> reach(N);
    vector<int> pred(N);
    vector<Edge*> pred_edge(N);
    vector<bool> visited(N);
    vector<vector<int>> children(N);

    SkarfQueue::MinIDQueue queue(N);
    int src = center->id;
    dist[src] = 0;
    queue.push({src, 0});

    while (!queue.empty()) {
        auto [vi, v_dist] = queue.pop();
        assert(v_dist == dist[vi]);
        if (vi != src) {
            assert(dist[vi] == dist[pred[vi]] + pred_edge[vi]->cost);
            children[pred[vi]].push_back(vi);
        }
        auto& neighbors = transposed ? nodes[vi]->back_adj : nodes[vi]->adj;
        for (Edge* e : neighbors) {
            int wi = e->other(vi);
            if (visited[wi]) continue;
            if (v_dist + e->cost < dist[wi]) {
                dist[wi] = v_dist + e->cost;
                pred[wi] = vi;
                pred_edge[wi] = e;
                if (queue.contains_id(wi))
                    queue.decrease_key({wi, dist[wi]});
                else
                    queue.push({wi, dist[wi]});
            }
        }
        visited[vi] = true;
    }

    // compute cut size for each step
    SegmentTree union_tree(num_steps + 1);
    auto& edges = graph->get_edges();
    for (int ei : skeleton_union) {
        int start = transposed? edges[ei]->to->id : edges[ei]->from->id;
        int left = dist[start] / step_size + 1;
        int right = (dist[start] + edges[ei]->cost) / step_size;
        union_tree.update(left, right + 1, 1);
    }
    vector<int> union_cut_sizes(num_steps + 1);
    for (int i = 0; i <= num_steps; ++i)
        union_cut_sizes[i] = union_tree.query(i, i + 1);

    // calculate maximum cut size of skeleton union from center node
    vector<pair<long long, long long>> intervals;
    for (int ei : skeleton_union) {
        int start = transposed? edges[ei]->to->id : edges[ei]->from->id;
        long long l = dist[start];
        long long r = l + edges[ei]->cost;
        intervals.emplace_back(l, r);
    }
    int maximum_union_cut_size = maximum_overlap(intervals);

    // compute edges in the skeleton from the center node
    vector<int> center_skeleton;
    vector<bool> visited_reach_computation(N);
    stack<int> s;
    s.emplace(src);
    while (!s.empty()) {
        int vi = s.top();
        if (!visited_reach_computation[vi]) {
            assert(reach[vi] == 0);
            for (int wi : children[vi]) {
                s.emplace(wi);
            }
            visited_reach_computation[vi] = true;
        } else {
            s.pop();
            // set the flags
            for (int wi : children[vi]) {
                Edge* e = pred_edge[wi];
                assert(reach[wi] + e->cost > 0);
                assert(dist[wi] == dist[vi] + e->cost);
                reach[vi] = max(reach[vi], e->cost + reach[wi]);
                // see also our skeleton definition for non-geometric-realisations
                if ((reach[wi] + e->cost) > dist[vi]) {
                    center_skeleton.push_back(e->id);
                }
            }
        }
    }

    // compute cut sizes of center skeleton
    SegmentTree center_tree(num_steps + 1);
    for (int ei : center_skeleton) {
        int start = transposed? edges[ei]->to->id : edges[ei]->from->id;
        int left = dist[start] / step_size + 1;
        int right = (dist[start] + edges[ei]->cost) / step_size;
        center_tree.update(left, right + 1, 1);
    }
    vector<int> center_cut_sizes(num_steps + 1);
    for (int i = 0; i <= num_steps; ++i)
        center_cut_sizes[i] = center_tree.query(i, i + 1);

    // calculate maximum cut size of spt from center node
    intervals.clear();
    for (int ei : center_skeleton) {
        int start = transposed? edges[ei]->to->id : edges[ei]->from->id;
        long long l = dist[start];
        long long r = l + edges[ei]->cost;
        intervals.emplace_back(l, r);
    }
    int maximum_center_cut_size = maximum_overlap(intervals);

    // export overlap values
    string folder = results_folder + "overlap/" + graph->region + "/" + to_string(partition_size) + "/" + to_string(cell_idx) + "/";
    validate_directory(folder);

    string prefix = folder;

    ofstream out(folder + "union_cuts.csv");
    out << "dist;count" << endl;
    for (int i = 0; i <= num_steps; ++i) {
        out << i * step_size << ";";
        out << union_cut_sizes[i] << endl;
    }
    out.close();

    out.open(folder + "center_cuts.csv");
    out << "dist;count" << endl;
    for (int i = 0; i <= num_steps; ++i) {
        out << i * step_size << ";";
        out << center_cut_sizes[i] << endl;
    }
    out.close();

    out.open(folder + "maximum_union_cut.txt");
    out << maximum_union_cut_size << endl;
    out.close();

    out.open(folder + "maximum_center_cut.txt");
    out << maximum_center_cut_size << endl;
    out.close();
}

void Skarf::precompute_cell(int cell_idx, int partition_size) {
    unordered_map<Node*, long long> boundary_dist_transposed;
    for (Node* node : graph->get_parts()[cell_idx]) boundary_dist_transposed[node] = 0;

    precompute_direction(cell_idx, true, partition_size, boundary_dist_transposed);
    Node* forward_center = min_element(boundary_dist_transposed.begin(), boundary_dist_transposed.end(),
                                       [](const auto& x, const auto& y) {
                                           return x.second < y.second;
                                       })
                               ->first;

    unordered_map<Node*, long long> boundary_dist;
    for (Node* node : graph->get_parts()[cell_idx]) boundary_dist[node] = 0;
    precompute_direction(cell_idx, false, partition_size, boundary_dist);
    Node* backward_center = min_element(boundary_dist.begin(), boundary_dist.end(),
                                        [](const auto& x, const auto& y) {
                                            return x.second < y.second;
                                        })
                                ->first;

    add_cell_overlap(cell_idx + partition_size, false, partition_size, forward_center);
    add_cell_overlap(cell_idx, true, partition_size, backward_center);
}

string round_to_precision(float value, int precision) {
    float p = pow(10, precision);
    string num = to_string(round(value * p) / p);
    int point_index = -1;
    for (int i = 0; i < int(num.size()); ++i)
        if (num[i] == '.') point_index = i;
    if (point_index != -1) num = num.substr(0, min(int(num.size()), point_index + 4));
    return num;
}

string get_memory_usage_in_gigabytes() {
    pid_t pid = getpid();
    string c = string("pmap ") + to_string(pid) + " | grep total";
    const char* cmd = c.c_str();
    array<char, 128> buffer;
    string result;
    unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    auto parsed = split(result, " ");
    string last = parsed[parsed.size() - 1];
    string kilobytes = last.substr(0, last.size() - 2);
    float gigabytes = atoi(kilobytes.c_str()) / 1e6;
    return round_to_precision(gigabytes, 3);
}

void Skarf::set_step_size() {
    int N = graph->get_n();
    auto& nodes = graph->get_nodes();

    vector<long long> dist(N, numeric_limits<long long>::max());
    vector<long long> reach(N);
    vector<int> pred(N);
    vector<Edge*> pred_edge(N);
    vector<bool> visited(N);
    vector<vector<int>> children(N);

    SkarfQueue::MinIDQueue queue(N);
    int src = 0;
    dist[src] = 0;
    queue.push({src, 0});
    long long largest_dist = 0;

    while (!queue.empty()) {
        auto [vi, v_dist] = queue.pop();
        assert(v_dist == dist[vi]);
        if (vi != src) {
            assert(dist[vi] == dist[pred[vi]] + pred_edge[vi]->cost);
            children[pred[vi]].push_back(vi);
        }
        largest_dist = max(largest_dist, v_dist);
        auto& neighbors = nodes[vi]->adj;
        for (Edge* e : neighbors) {
            int wi = e->other(vi);
            if (visited[wi]) continue;
            if (v_dist + e->cost < dist[wi]) {
                dist[wi] = v_dist + e->cost;
                pred[wi] = vi;
                pred_edge[wi] = e;
                if (queue.contains_id(wi))
                    queue.decrease_key({wi, dist[wi]});
                else
                    queue.push({wi, dist[wi]});
            }
        }
        visited[vi] = true;
    }

    long long max_dist = largest_dist * 2;
    step_size = max_dist / num_steps;
}

void Skarf::precompute(int start, int end) {
    auto precomputation_start_time = chrono::high_resolution_clock::now();

    set_step_size();

    synced_stream sync_out;
    thread_pool pool(this->num_threads);

    string finished_path = graph->region_folder + "/intermediate/finished_" + to_string(partition_size) + "_" + get_graph()->get_current_time_point_as_string() + ".csv";
    ofstream finished(finished_path);

    for (int i = start; i < end; i++) {
        pool.push_task(precompute_cell, i, partition_size);
    }

    while (pool.get_tasks_total() > 0) {
        sync_out.println(pool.get_tasks_total(),
                         " tasks total, ",
                         pool.get_tasks_running(),
                         " tasks running, ",
                         pool.get_tasks_queued(),
                         " tasks queued, ",
                         get_memory_usage_in_gigabytes(),
                         " GB memory usage");
        if (pool.get_tasks_total() < partition_size) {
            auto currentTime = chrono::high_resolution_clock::now();
            size_t minutesPassed = chrono::duration_cast<chrono::minutes>(currentTime - precomputation_start_time).count();
            long long hoursPassed = (long long)(minutesPassed / 60);
            long long minutesPassedMod60 = minutesPassed % 60;
            long long minutesLeft = (long long)(minutesPassed * pool.get_tasks_total() / (partition_size - pool.get_tasks_total()));
            long long hoursLeft = (long long)(minutesLeft / 60);
            long long minutesLeftMod60 = minutesLeft % 60;
            sync_out.println("Finished ",
                             (partition_size - pool.get_tasks_total()),
                             " cells in ",
                             hoursPassed,
                             ":",
                             minutesPassedMod60,
                             " hours. Approximately ",
                             hoursLeft,
                             ":",
                             minutesLeftMod60,
                             " hours left.");
        }
        this_thread::sleep_for(std::chrono::milliseconds(10000));
    }
    pool.wait_for_tasks();
    precomputed = true;

    merge_flags("arcflags");
    merge_flags("skarf");

    auto finished_cell_time = chrono::high_resolution_clock::now();
    size_t hours = chrono::duration_cast<chrono::hours>(finished_cell_time - precomputation_start_time).count();
    size_t minutes = chrono::duration_cast<chrono::minutes>(finished_cell_time - precomputation_start_time).count() % 60;
    size_t seconds = chrono::duration_cast<chrono::seconds>(finished_cell_time - precomputation_start_time).count() % 60;
    finished << "finished flag precomputation in " << hours << ":" << minutes << ":" << seconds << " hours" << endl;
    cout << "finished flag precomputation in " << hours << ":" << minutes << ":" << seconds << " hours" << endl;
    finished.close();
}

void Skarf::compress(unordered_map<Edge*, BitLabel>& labels) {
    preprocessing_edge_to_key = vector<size_t>(graph->get_m());
    unordered_map<size_t, BitLabel> temp_labels;
    hash<BitLabel> hash_f;
    for (auto& [edge, label] : labels) {
        size_t key = hash_f(label);
        if (temp_labels.find(key) != temp_labels.end()) {
            if (label != temp_labels[key])
                cout << "Hash collision!" << endl;
        } else
            temp_labels[key] = label;
        preprocessing_edge_to_key[edge->id] = key;
    }
    temp_labels[0].resize(2 * partition_size, 0);

    preprocessing_key_to_label = vector<BitLabel>(temp_labels.size());
    int i = 0;
    unordered_map<size_t, int> id_map;
    for (auto& [key, label] : temp_labels) {
        id_map[key] = i;
        preprocessing_key_to_label[i] = label;
        ++i;
    }
    for (int edge = 0; edge < preprocessing_edge_to_key.size(); ++edge) {
        auto key = preprocessing_edge_to_key[edge];
        preprocessing_edge_to_key[edge] = id_map[key];
    }
}

void Skarf::merge_flags(string folder) {
    unordered_map<Edge*, BitLabel> labels;
    for (int i = 0; i < 2 * partition_size; i++) {
        if (i % 10 == 0) cout << "merge cell " << i + 1 << " of " << 2 * partition_size << endl;
        string dir = graph->region_folder + "/intermediate/" + folder;
        validate_directory(dir);
        ifstream file(dir + "/flags_partitionSize" + to_string(partition_size) + "_cell" + to_string(i));
        string line;
        while (getline(file, line)) {
            long long edge_id = stoll(line);
            Edge* edge = get_graph()->get_edges()[edge_id];
            if (labels.find(edge) == labels.end()) labels[edge].resize(2 * partition_size, 0);
            labels[edge][i] = 1;
        }
    }
    cout << "Start compressing" << endl;
    compress(labels);
    cout << "Start exporting" << endl;
    export_flags(folder);
    cout << "Finished exporting" << endl;
    precomputed = true;
}

void Skarf::import_all_flags(bool import_arc_flags, bool import_skarf_flags) {
    if (import_arc_flags) {
        cout << "import arc flags" << endl;
        import_flags(0);
        cout << "finished importing arc flags" << endl;
    }
    if (import_skarf_flags) {
        cout << "import skarf flags" << endl;
        import_flags(1);
        cout << "finished importing skarf flags" << endl;
    }
}

void Skarf::import_flags(int flagType) {
    auto& key_to_label = flagType ? skarf_key_to_label : arcflags_key_to_label;
    vector<size_t>& edge_to_key = flagType ? skarf_edge_to_key : arcflags_edge_to_key;
    edge_to_key = vector<size_t>(graph->get_m());

    string flag_type = flagType ? "skarf" : "arcflags";
    string dir = graph->region_folder + "/flags/";
    validate_directory(dir);
    string flags_path = dir + flag_type + "_" + to_string(partition_size) + "_hashToFlag.bin";
    string edges_path = dir + flag_type + "_" + to_string(partition_size) + "_edgeToHash.csv";

    ifstream file(edges_path);
    string line;
    getline(file, line);  // skip first line
    vector<string> lines;
    size_t max_key = 0;
    while (getline(file, line)) {
        vector<string> csv = split(line, ";", false);
        long long edge_id = stoll(csv[0]);
        size_t key = stoul(csv[1]);
        if (key > max_key) max_key = key;
        edge_to_key[edge_id] = key;
    }

    key_to_label.clear();
    key_to_label.resize(max_key + 1);

    ifstream input(flags_path, ios::binary);
    vector<unsigned char> buffer(istreambuf_iterator<char>(input), {});
    int row_byte_size = sizeof(size_t) + ceil(2 * partition_size / 8.0);

    bitset<64> key_bits;
    BitLabel flag(2 * partition_size);
    size_t key;
    for (int i = 0; i < buffer.size(); i++) {
        int pos = i % row_byte_size;
        if (pos == 0) {
            key_bits = bitset<64>{0};
            flag = BitLabel(2 * partition_size, 0);
        }
        if (pos < 8) {  // read key
            bitset<64> mask{buffer[i]};
            key_bits |= mask;
            if (pos == 7) {
                key = key_bits.to_ulong();
            } else
                key_bits <<= 8;
        } else {  // read flag
            flag |= BitLabel(2 * partition_size, buffer[i]);
            if (pos < row_byte_size - 1) {
                flag <<= 8;
            } else
                key_to_label[key] = flag;
        }
    }
    precomputed = true;
}

void Skarf::export_flags(string flag_type) {
    string dir = graph->region_folder + "/flags/";
    validate_directory(dir);
    ofstream csv_file(dir + flag_type + "_" + to_string(partition_size) + "_edgeToHash.csv");
    csv_file << "edge_id;key\n";
    for (int i = 0; i < preprocessing_edge_to_key.size(); i++) {
        csv_file << i << ";" << preprocessing_edge_to_key[i] << "\n";
    }
    csv_file.close();

    vector<byte> bytes;
    for (int i = 0; i < preprocessing_key_to_label.size(); ++i) {
        size_t key = i;
        auto label = preprocessing_key_to_label[i];

        BitLabel flag = label;
        bitset<64> key_bits(key_bits);
        vector<byte> key_bytes;
        for (int i = 0; i < sizeof(size_t); i++) {
            byte nextByte{static_cast<byte>(key & ((1 << 8) - 1))};
            key_bytes.push_back(nextByte);
            key >>= 8;
        }
        for (int i = 0; i < key_bytes.size(); i++)
            bytes.push_back(key_bytes[key_bytes.size() - 1 - i]);

        vector<byte> flag_bytes;
        for (int i = 0; i < 2 * partition_size; i += 8) {
            BitLabel copy = flag;
            copy &= BitLabel(2 * partition_size, 255);
            byte nextByte{static_cast<byte>(copy.to_ulong())};
            flag_bytes.push_back(nextByte);
            flag >>= 8;
        }
        for (int i = 0; i < flag_bytes.size(); i++)
            bytes.push_back(flag_bytes[flag_bytes.size() - 1 - i]);
    }
    ofstream outfile(dir + flag_type + "_" + to_string(partition_size) + "_hashToFlag.bin", ios::out | ios::binary);
    outfile.write((const char*)&bytes[0], bytes.size());

    backup_flags(flag_type);
}

void Skarf::backup_flags(string flag_type) {
    string dir = graph->region_folder;
    validate_directory(dir + "/flags/");
    validate_directory(dir + "/flags_backup/");
    string time_stamp = graph->get_current_time_point_as_string();
    vector<pair<string, string>> types = {{"edgeToHash", ".csv"}, {"hashToFlag", ".bin"}};
    for (auto [type, suffix] : types) {
        string src = dir + "/flags/" + flag_type + "_" + to_string(partition_size) + "_" + type + suffix;
        string dest = dir + "/flags_backup/" + flag_type + "_" + to_string(partition_size) + "_" + type + "_" + time_stamp + suffix;
        system(("cp " + src + " \"" + dest + "\"").c_str());
    }
}

bool Skarf::is_set_skarf(int edge, int startCell, int destination_cell, int flagType) {
    vector<size_t>& edge_to_key = flagType ? skarf_edge_to_key : arcflags_edge_to_key;
    if (!edge_to_key[edge]) return false;
    auto& edge_label_hash = edge_to_key[edge];
    auto& edge_label = flagType ? skarf_key_to_label[edge_label_hash] : arcflags_key_to_label[edge_label_hash];
    if (flagType == 0) return edge_label[destination_cell];
    return edge_label[destination_cell] || edge_label[startCell + partition_size];
}

bool Skarf::is_set_skarf_plus(int edge, int startCell, int destination_cell) {
    auto& skarf_label = skarf_key_to_label[skarf_edge_to_key[edge]];
    auto& arcflags_label = arcflags_key_to_label[arcflags_edge_to_key[edge]];
    return arcflags_label[destination_cell] && (skarf_label[destination_cell] || skarf_label[startCell + partition_size]);
}

bool Skarf::is_set_arcflags(int edge, int destination_cell) {
    auto& arcflags_label = arcflags_key_to_label[arcflags_edge_to_key[edge]];
    return arcflags_label[destination_cell];
}

bool Skarf::is_set_skarf_plus_forward(int edge, int startCell, int destination_cell) {
    auto& skarf_label = skarf_key_to_label[skarf_edge_to_key[edge]];
    auto& arcflags_label = arcflags_key_to_label[arcflags_edge_to_key[edge]];
    return arcflags_label[destination_cell] && skarf_label[startCell + partition_size];
}

bool Skarf::is_set_skarf_plus_backward(int edge, int startCell, int destination_cell) {
    auto& skarf_label = skarf_key_to_label[skarf_edge_to_key[edge]];
    auto& arcflags_label = arcflags_key_to_label[arcflags_edge_to_key[edge]];
    return arcflags_label[startCell + partition_size] && skarf_label[destination_cell];
}

bool Skarf::is_set_arcflags_forward(int edge, int startCell, int destination_cell) {
    auto& arcflags_label = arcflags_key_to_label[arcflags_edge_to_key[edge]];
    return arcflags_label[destination_cell];
}

bool Skarf::is_set_arcflags_backward(int edge, int startCell, int destination_cell) {
    auto& arcflags_label = arcflags_key_to_label[arcflags_edge_to_key[edge]];
    return arcflags_label[startCell + partition_size];
}

list<Edge*> Skarf::calculate_route_with(int source_id, int target_id, bool arcflags, bool skarf) {
    int N = get_graph()->get_n();
    auto& nodes = get_graph()->get_nodes();
    Node* source = nodes[source_id];
    Node* target = nodes[target_id];
    int source_cell = source->cell_idx;
    int target_cell = target->cell_idx;

    unordered_map<Node*, long long> dist;
    unordered_map<Node*, Node*> pred;
    unordered_map<Node*, Edge*> pred_edge;
    unordered_map<Node*, bool> visited;

    dist[source] = 0;
    forward_queue.push({source_id, 0});

    while (!forward_queue.empty()) {
        auto [vi, v_dist] = forward_queue.pop();
        if (vi == target_id) break;
        Node* v = nodes[vi];
        settle_node();
        for (Edge* e : v->adj) {
            visit_edge();  // count edges that are looked at
            Node* w = e->to;
            if (visited[w]) continue;
            if (arcflags) {
                if (skarf) {
                    if (!is_set_skarf_plus(e->id, source_cell, target_cell)) continue;
                } else {
                    if (!is_set_arcflags(e->id, target_cell)) continue;
                }
            }
            if (dist[w] == 0 || v_dist + e->cost < dist[w]) {
                long long w_dist = v_dist + e->cost;
                dist[w] = w_dist;
                pred[w] = v;
                pred_edge[w] = e;
                if (forward_queue.contains_id(w->id))
                    forward_queue.decrease_key({w->id, w_dist});
                else
                    forward_queue.push({w->id, w_dist});
            }
        }
        visited[v] = true;
    }

    list<Edge*> edges;
    Node* curr = target;
    while (pred_edge[curr]) {
        edges.push_front(pred_edge[curr]);
        curr = pred[curr];
    }
    return edges;
}

list<Edge*> Skarf::calculate_route_bidirectional_with(int source_id, int target_id, bool arcflags, bool skarf) {
    int N = get_graph()->get_n();
    auto& nodes = get_graph()->get_nodes();
    Node* source = nodes[source_id];
    Node* target = nodes[target_id];
    int source_cell = source->cell_idx;
    int target_cell = target->cell_idx;

    unordered_map<Node*, long long> dist_forward;
    unordered_map<Node*, Node*> pred_forward;
    unordered_map<Node*, Edge*> pred_edge_forward;
    unordered_map<Node*, bool> visited_forward;

    unordered_map<Node*, long long> dist_backward;
    unordered_map<Node*, Node*> pred_backward;
    unordered_map<Node*, Edge*> pred_edge_backward;
    unordered_map<Node*, bool> visited_backward;

    dist_forward[source] = 0;
    dist_backward[target] = 0;
    forward_queue.push({source_id, 0});
    backward_queue.push({target_id, 0});
    Node* via_node;
    long long min_cost = numeric_limits<long long>::max();

    while (!forward_queue.empty() || !backward_queue.empty()) {
        bool search_forward;
        if (forward_queue.empty())
            search_forward = false;
        else if (backward_queue.empty())
            search_forward = true;
        else
            search_forward = forward_queue.peek().key <= backward_queue.peek().key;

        if (search_forward) {
            auto [vi, v_dist] = forward_queue.pop();
            if (vi == target_id) break;
            Node* v = nodes[vi];
            settle_node();
            for (Edge* e : v->adj) {
                visit_edge();
                Node* w = e->to;
                if (visited_forward[w]) continue;
                if (arcflags) {
                    if (skarf) {
                        if (!is_set_skarf_plus_forward(e->id, source_cell, target_cell)) continue;
                    } else {
                        if (!is_set_arcflags_forward(e->id, source_cell, target_cell)) continue;
                    }
                }
                if (dist_forward[w] == 0 || v_dist + e->cost < dist_forward[w]) {
                    long long w_dist = v_dist + e->cost;
                    dist_forward[w] = w_dist;
                    pred_forward[w] = v;
                    pred_edge_forward[w] = e;
                    if (forward_queue.contains_id(w->id))
                        forward_queue.decrease_key({w->id, w_dist});
                    else
                        forward_queue.push({w->id, w_dist});
                    if (visited_backward[w]) {
                        long long via_dist = w_dist + dist_backward[w];
                        if (via_dist < min_cost) {
                            min_cost = via_dist;
                            via_node = w;
                        }
                    }
                }
            }
            visited_forward[v] = true;
            if (visited_backward[v] == true) {
                break;
            }
        } else {
            auto [vi, v_dist] = backward_queue.pop();
            if (vi == source_id) break;
            Node* v = nodes[vi];
            settle_node();
            for (Edge* e : v->back_adj) {
                visit_edge();
                Node* w = e->from;
                if (visited_backward[w]) continue;
                if (arcflags) {
                    if (skarf) {
                        if (!is_set_skarf_plus_backward(e->id, source_cell, target_cell)) continue;
                    } else {
                        if (!is_set_arcflags_backward(e->id, source_cell, target_cell)) continue;
                    }
                }
                if (dist_backward[w] == 0 || v_dist + e->cost < dist_backward[w]) {
                    long long w_dist = v_dist + e->cost;
                    dist_backward[w] = w_dist;
                    pred_backward[w] = v;
                    pred_edge_backward[w] = e;
                    if (backward_queue.contains_id(w->id))
                        backward_queue.decrease_key({w->id, w_dist});
                    else
                        backward_queue.push({w->id, w_dist});
                    if (visited_forward[w]) {
                        long long via_dist = dist_forward[w] + w_dist;
                        if (via_dist < min_cost) {
                            min_cost = via_dist;
                            via_node = w;
                        }
                    }
                }
            }
            visited_backward[v] = true;
            if (visited_forward[v] == true) {
                break;
            }
        }
    }
    list<Edge*> edges;
    Node* curr = via_node;
    while (pred_forward[curr] != 0) {
        edges.push_front(pred_edge_forward[curr]);
        curr = pred_forward[curr];
    }
    curr = via_node;
    while (pred_backward[curr] != 0) {
        edges.push_back(pred_edge_backward[curr]);
        curr = pred_backward[curr];
    }
    return edges;
}

void Skarf::set_graph(Graph* g) {
    graph = g;
    forward_queue = SkarfQueue::MinIDQueue(g->get_n());
    backward_queue = SkarfQueue::MinIDQueue(g->get_n());
}

void Skarf::reset() {
    visited_edges = 0;
    settled_nodes = 0;
    forward_queue.clear();
    backward_queue.clear();
}

Graph* Skarf::graph;
SkarfQueue::MinIDQueue Skarf::forward_queue;
SkarfQueue::MinIDQueue Skarf::backward_queue;
