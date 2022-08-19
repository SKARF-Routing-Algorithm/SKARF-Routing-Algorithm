#include <cxxopts.h>
#include <memory.h>

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "boost/filesystem.hpp"
#include "experiments.h"
#include "graph.h"
#include "skarf.h"
#include "util.h"

using namespace std;
bool _compare_with_dijkstra;

string load_config() {
    ifstream in("../data_directory.txt");
    string data_dir;
    in >> data_dir;
    if (data_dir.back() != '/') data_dir += '/';
    return "../" + data_dir;
}

void init_skarf(Graph& g, string dir, unsigned partition_size, unsigned runs, unsigned threads, string data_dir, bool precompute = false) {
    string results_folder = data_dir.substr(0, data_dir.size() - 5) + "results/";
    Experiments engine = Experiments(&g);
    engine.results_folder = results_folder;

    Skarf::get_instance()->results_folder = results_folder;
    Skarf::get_instance()->set_graph(&g);
    Skarf::get_instance()->set_partition_size(partition_size);
    Skarf::get_instance()->set_threads(threads);
    g.partition(partition_size);
    for (auto type : {"arcflags", "skarf"}) {
        ifstream csv(dir + type + "_" + to_string(partition_size) + "_edgeToHash.csv");
        ifstream bin(dir + type + "_" + to_string(partition_size) + "_hashToFlag.bin");
        if (!csv || !bin) precompute = true;
    }
    if (precompute) {
        cout << "Start SKARF computation." << endl;
        Skarf::get_instance()->precompute(0, partition_size);
        cout << "Finished SKARF computation." << endl;
    }
    if (runs > 0) {
        cout << "Start comparisons" << endl;
        Skarf::get_instance()->import_all_flags(true, true);
        engine.run_comparisons(runs, partition_size);
        cout << "Finished comparisons" << endl;
    }
}

int parse(int argc, const char* argv[]) {
    try {
        unique_ptr<cxxopts::Options> allocated(new cxxopts::Options(argv[0], " - Routing command line options"));
        auto& options = *allocated;
        options
            .positional_help("[optional args]")
            .show_positional_help();
        options.set_width(90).set_tab_expansion().allow_unrecognised_options().add_options()
            ("h,help", "Print help")
            ("g,graph", "Input graph (e.g. germany, ...)", cxxopts::value<string>(), "name")
            ("r,runs", "Number of test queries", cxxopts::value<unsigned>()->default_value("0"), "number")
            ("p,partition", "Partition size", cxxopts::value<unsigned>()->default_value("128"), "size")
            ("t,threads", "Number of parallel threads", cxxopts::value<unsigned>()->default_value("4"), "number")
            ("v,validate", "Validate query results against dijkstra")
            ("f,force", "Force flag precomputation (overwrite existing flags)");
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            cout << options.help({"", "CHASE specific"}) << endl;
            return true;
        }

        Graph g;
        unsigned num_threads = result["threads"].as<unsigned>();
        unsigned partition_size = result["partition"].as<unsigned>();
        unsigned runs = result["runs"].as<unsigned>();
        bool precompute = result.count("force");
        string region = result["graph"].as<string>();
        g.region = region;
        string data_dir = load_config();

        // load graph
        if (result.count("graph")) {
            validate_directory(data_dir + region);
            g.load_network(data_dir + region);
        } else {
            cout << "Please provide a graph file! (-g[--graph] name)" << endl;
            return false;
        }

        if (result.count("validate")) {
            _compare_with_dijkstra = true;
        }
        init_skarf(g, data_dir + region + "/flags/", partition_size, runs, num_threads, data_dir, precompute);

    } catch (const cxxopts::exceptions::exception& e) {
        cout << "error parsing options: " << e.what() << endl;
        return false;
    }
    return true;
}

int main(int argc, const char* argv[]) {
    if (!parse(argc, argv)) {
        return 1;
    }
    return 0;
}