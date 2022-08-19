#include "util.h"

#include <cmath>

using namespace std;

bool contains(string& delimiters, char c) {
    bool result = false;
    for (char ch : delimiters) result |= ch == c;
    return result;
}

vector<string> split(string& line, string delimiters, bool empty_allowed) {
    vector<string> tokens;
    size_t start = 0, end = 0;
    while (start < line.size()) {
        while (end < line.size() && !contains(delimiters, line[end])) ++end;
        if (start < end || empty_allowed)
            tokens.push_back(line.substr(start, end - start));
        start = ++end;
    }
    return tokens;
}

void validate_directory(string dir) {
    if (!filesystem::is_directory(dir))
        system(("mkdir -p " + dir).c_str());
}

ifstream read_file(string path){
    auto parts = split(path, "/");
    if (parts.size() > 1) {
        string dir = path.substr(0, path.length() - ("/" + parts.back()).size());
        validate_directory(dir);
    }
    return ifstream(path);
}
