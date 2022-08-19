#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>

using namespace std;

vector<string> split(string& line, string delim, bool empty_allowed = true);
void validate_directory(string path);
ifstream read_file(string path);