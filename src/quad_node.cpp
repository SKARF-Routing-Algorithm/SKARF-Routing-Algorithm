#include "quad_node.h"

#include <algorithm>
#include <queue>

#include "point.h"

using namespace std;

QuadNode::QuadNode(Point* topL, Point* botR) {
    node = NULL;
    filled = false;
    top_left_tree = NULL;
    top_right_tree = NULL;
    bot_left_tree = NULL;
    bot_right_tree = NULL;
    topLeft = topL;
    bot_right = botR;
}

void QuadNode::insert(struct Node* node) {
    if (node == NULL)
        return;

    if (!contains(Point(node->lon, node->lat)))
        return;

    if (!filled) {
        this->node = node;
        filled = true;
        return;
    }

    if (filled && this->node != 0) {
        struct Node* temp = this->node;
        this->node = NULL;
        insert(temp);
    }

    if ((topLeft->get_lon() + bot_right->get_lon()) / 2 >= node->lon) {
        // Indicates topLeftTree
        if ((topLeft->get_lat() + bot_right->get_lat()) / 2 <= node->lat) {
            if (top_left_tree == NULL) {
                top_left_tree = new QuadNode(
                    topLeft,
                    new Point((topLeft->get_lon() + bot_right->get_lon()) / 2,
                              (topLeft->get_lat() + bot_right->get_lat()) / 2));
            }
            top_left_tree->insert(node);
        }
        // Indicates botLeftTree
        else {
            if (bot_left_tree == NULL) {
                bot_left_tree = new QuadNode(
                    new Point(topLeft->get_lon(),
                              (topLeft->get_lat() + bot_right->get_lat()) / 2),
                    new Point((topLeft->get_lon() + bot_right->get_lon()) / 2,
                              bot_right->get_lat()));
            }
            bot_left_tree->insert(node);
        }
    } else {
        // Indicates topRightTree
        if ((topLeft->get_lat() + bot_right->get_lat()) / 2 <= node->lat) {
            if (top_right_tree == NULL) {
                top_right_tree = new QuadNode(
                    new Point((topLeft->get_lon() + bot_right->get_lon()) / 2,
                              topLeft->get_lat()),
                    new Point(bot_right->get_lon(),
                              (topLeft->get_lat() + bot_right->get_lat()) / 2));
            }
            top_right_tree->insert(node);
        }
        // Indicates botRightTree
        else {
            if (bot_right_tree == NULL) {
                bot_right_tree = new QuadNode(
                    new Point((topLeft->get_lon() + bot_right->get_lon()) / 2,
                              (topLeft->get_lat() + bot_right->get_lat()) / 2),
                    bot_right);
            }
            bot_right_tree->insert(node);
        }
    }
}

Node* QuadNode::search(Point* p) {
    cout << "search in [(" << topLeft->get_lon() << "," << topLeft->get_lat() << "), (" << bot_right->get_lon() << ", " << bot_right->get_lat() << ")]" << endl;
    if (!contains(*p))
        return NULL;

    if (this->node != 0)
        return this->node;

    if ((topLeft->get_lon() + bot_right->get_lon()) / 2 >= p->get_lon()) {
        // Indicates topLeftTree
        if ((topLeft->get_lat() + bot_right->get_lat()) / 2 <= p->get_lat()) {
            if (top_left_tree == NULL)
                return NULL;
            top_left_tree->search(p);
        }
        // Indicates botLeftTree
        else {
            if (bot_left_tree == NULL)
                return NULL;
            bot_left_tree->search(p);
        }
    } else {
        // Indicates topRightTree
        if ((topLeft->get_lat() + bot_right->get_lat()) / 2 <= p->get_lat()) {
            if (top_right_tree == NULL)
                return NULL;
            top_right_tree->search(p);
        }
        // Indicates botRightTree
        else {
            if (bot_right_tree == NULL)
                return NULL;
            bot_right_tree->search(p);
        }
    }
    return NULL;
};

double QuadNode::euclidean_distance(double lat1, double lon1, double lat2, double lon2) {
    double delta_lat = abs(lat1 - lat2);
    double delta_lon = abs(lon1 - lon2);
    double x = sqrt(pow(delta_lat, 2) + pow(delta_lon, 2));
    return x;
}

struct Result QuadNode::nearest(Point* point, Result best) {
    double x1 = this->topLeft->get_lon();
    double y1 = this->topLeft->get_lat();
    double x2 = this->bot_right->get_lon();
    double y2 = this->bot_right->get_lat();

    double x = point->get_lon();
    double y = point->get_lat();
    if (x < x1 - best.dist || x > x2 + best.dist || y > y1 + best.dist || y < y2 - best.dist)
        return best;
    auto n = this->node;
    if (n) {
        double d = euclidean_distance(y, x, n->lat, n->lon);
        if (d < best.dist) {
            best.dist = d;
            best.node = n;
        }
    }
    if (top_left_tree) best = top_left_tree->nearest(point, best);
    if (top_right_tree) best = top_right_tree->nearest(point, best);
    if (bot_left_tree) best = bot_left_tree->nearest(point, best);
    if (bot_right_tree) best = bot_right_tree->nearest(point, best);
    return best;
}

QuadNode* QuadNode::get_child(int num) {
    switch (num) {
        case 0:
            return top_left_tree;
        case 1:
            return top_right_tree;
        case 2:
            return bot_left_tree;
        case 3:
            return bot_right_tree;
    }
    return NULL;
}

bool QuadNode::contains(Point p) {
    return (p.get_lon() >= topLeft->get_lon() &&
            p.get_lat() <= topLeft->get_lat() &&
            p.get_lon() <= bot_right->get_lon() &&
            p.get_lat() >= bot_right->get_lat());
}