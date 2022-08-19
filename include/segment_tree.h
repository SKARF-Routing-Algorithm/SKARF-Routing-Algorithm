#include<vector>
using namespace std;

// Segment tree with right exclusive intervals
struct SegmentTree {

private:
	using T = int; // data type
	const T NEUTRAL = 0;
	using U = int; // lazy type
	const U NO_UPD = 0;
	vector<T> seg; // aggregates
	vector<U> lazy;
	int m_n;
	T combine(T a, T b) { return a + b; }
	T apply(int v, const U& upd) {
		seg[v] += upd;
		lazy[v] += upd;
		return seg[v];
	}
	void push(int idx, int l, int r) {
		int m = (r + l) / 2;
		if (m - l) apply(2 * idx, lazy[idx]);
		if (r - m) apply(2 * idx + 1, lazy[idx]);
		lazy[idx] = NO_UPD;
	}
	T query(int idx, int l, int r, int ql, int qr) {
		if (qr <= l || r <= ql) return NEUTRAL;
		if (ql <= l && r <= qr) return seg[idx];
		push(idx, l, r);
		return combine(query(2 * idx, l, (l + r) / 2, ql, qr), query(2 * idx + 1, (l + r) / 2, r, ql, qr));
	}
	T range_update(int idx, int l, int r, int ql, int qr, U upd) {
		if (qr <= l || r <= ql) return seg[idx];
		if (ql <= l && r <= qr) return apply(idx, upd);
		push(idx, l, r);
		return seg[idx] = combine(range_update(2 * idx, l, (l + r) / 2, ql, qr, upd), range_update(2 * idx + 1, (l + r) / 2, r, ql, qr, upd));
	}

public:
	SegmentTree(int n) : seg(4 * n, NEUTRAL), lazy(4 * n, NO_UPD), m_n(n) {};
	T query(int ql, int qr) { return query(1, 0, m_n, ql, qr); }
	void update(int ql, int qr, U upd) { range_update(1, 0, m_n, ql, qr, upd); }
};