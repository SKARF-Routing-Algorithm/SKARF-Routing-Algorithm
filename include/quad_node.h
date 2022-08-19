#pragma once
#include <cmath>
#include <iostream>
#include <vector>

#include "point.h"
using namespace std;

struct Edge;

struct Node {
    Node(int _id, double _lat, double _lon) : id{_id}, lat{_lat}, lon{_lon} {}
    int degree;
    int cell_idx;
    int id;
    bool boundaryNode;
    double lat, lon;
    string s_lat, s_lon;
    std::vector<Edge*> adj;
    std::vector<Edge*> back_adj;
};

struct Edge {
    Edge(Node* _from, Node* _to) : from{_from}, to{_to} {
    }
    int id;
    struct Node* from;
    struct Node* to;
    struct Edge* next;
    long long cost;
    double km, kmh;
    Node* other(Node* n) {
        return n->id == from->id ? to : from;
    }
    int other(int v) {
        return v == from->id ? to->id : from->id;
    }
};

struct Result {
    Result(Node* _node, double _dist) {
        node = _node;
        dist = _dist;
    }
    struct Node* node;
    double dist;
};

class QuadNode {
    // hold details of the boundary of this node
    Point* topLeft;
    Point* bot_right;

    // contains details of node
    struct Node* node;
    bool filled;

    // children of this tree
    QuadNode* top_left_tree;
    QuadNode* top_right_tree;
    QuadNode* bot_left_tree;
    QuadNode* bot_right_tree;

   public:
    QuadNode();
    QuadNode(Point* topL, Point* botR);
    void insert(struct Node*);
    struct Node* search(Point* point);
    bool contains(Point point);
    Result nearest(Point* point, Result best);
    double euclidean_distance(double lat1, double lon1, double lat2, double lon2);
    QuadNode* get_child(int num);
    Point* get_top_left() { return topLeft; };
    Point* get_bot_right() { return bot_right; }
};