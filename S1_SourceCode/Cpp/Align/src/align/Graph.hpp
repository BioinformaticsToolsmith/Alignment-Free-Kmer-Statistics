#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <algorithm>
#include <iostream>
#include <vector>
#include <tuple>
using namespace std;

// int transform(char c)
// {
// 	int ret = 0;
// 	int r;
// 	switch (c) {
// 	case 'a':
// 	case 'A':
// 		ret = 0;
// 		break;
// 	case 'c':
// 	case 'C':
// 		ret = 1;
// 		break;
// 	case 'g':
// 	case 'G':
// 		ret = 2;
// 		break;
// 	case 't':
// 	case 'T':
// 	case 'u':
// 	case 'U':
// 		ret = 3;
// 		break;
// 	case 'r':
// 	case 'R': // a or g
// 		if (rand() % 2) {
// 			ret = transform('A');
// 		} else {
// 			ret = transform('G');
// 		}
// 		break;
// 	case 'y':
// 	case 'Y': // c, t, or u
// 		r = rand() % 2;
// 		if (r == 0) {
// 			ret = transform('C');
// 		} else {
// 			ret = transform('T');
// 		}
// 		break;
// 	case 'k':
// 	case 'K': // g, t, or u
// 		r = rand() % 2;
// 		if (r == 0) {
// 			ret = transform('G');
// 		} else {
// 			ret = transform('T');
// 		}
// 		break;
// 	case 'm':
// 	case 'M':
// 		r = rand() % 2;
// 		if (r == 0) {
// 			ret = transform('A');
// 		} else {
// 			ret = transform('C');
// 		}
// 		break;
// 	case 's':
// 	case 'S':
// 		r = rand() % 2;
// 		if (r == 0) {
// 			ret = transform('C');
// 		} else {
// 			ret = transform('G');
// 		}
// 		break;
// 	case 'w':
// 	case 'W':
// 		r = rand() % 2;
// 		if (r == 0) {
// 			ret = transform('A');
// 		} else {
// 			ret = transform('T');
// 		}
// 		break;
// 	case 'b':
// 	case 'B':
// 		r = rand() % 3;
// 		if (r == 0) {
// 			ret = transform('C');
// 		} else if (r == 1) {
// 			ret = transform('G');
// 		} else {
// 			ret = transform('T');
// 		}
// 		break;
// 	case 'd':
// 	case 'D':
// 		r = rand() % 3;
// 		if (r == 0) {
// 			ret = transform('A');
// 		} else if (r == 1) {
// 			ret = transform('G');
// 		} else {
// 			ret = transform('T');
// 		}
// 		break;
// 	case 'h':
// 	case 'H':
// 		r = rand() % 3;
// 		if (r == 0) {
// 			ret = transform('A');
// 		} else if (r == 1) {
// 			ret = transform('C');
// 		} else {
// 			ret = transform('T');
// 		}
// 		break;
// 	case 'v':
// 	case 'V':
// 		r = rand() % 3;
// 		if (r == 0) {
// 			ret = transform('A');
// 		} else if (r == 1) {
// 			ret = transform('C');
// 		} else {
// 			ret = transform('G');
// 		}
// 		break;
// 	case 'n':
// 	case 'N':
// 		ret = rand() % 4;
// 		break;
// 	default:
// 		throw "bad nucl";
// 	}
// 	return ret;
// }
inline int transform(int c) { return c; };
class Graph {
public:
	Graph(string str1, string str2, const double matrix[4][4], double sigma_0, double epsilon_0, bool local_align=true)
		: s1(str1), s2(str2), sigma(sigma_0), epsilon(epsilon_0), l1(str1.length()+1), l2(str2.length()+1), local(local_align) { init(matrix); };
	Graph(string str1, string str2, double sigma_0, double epsilon_0, bool local_align=true) : s1(str1), s2(str2), sigma(sigma_0), epsilon(epsilon_0), l1(str1.length()+1), l2(str2.length()+1), local(local_align) {};
	~Graph() { delete[] graph; delete[] backtrack; };
	vector<string> get_alignment(std::pair<int,int> &begin, std::pair<int,int> &end) const;
	double get_total_score() const {
		return graph[hash(l1-1, l2-1, 1)];
	}
	double get_identity(vector<string> alignment) const;
	void init(const double mat[][4]);
private:
	inline int hash(int x, int y, int z) const { return z + 3 * (y + l2 * x); };

	double get_score(char c1, char c2) const;
        void fill_graph_at(int i, int j);
        double *graph;
	tuple<int,int,int> *backtrack;
	double score[4][4];
	const string s1,s2;
	const double sigma, epsilon;
	const int l1,l2;
	const bool local;
        tuple<int,int,int> local_end;
};

double Graph::get_identity(vector<string> alignment) const {
	int count = 0;
	for (int i = 0; i < alignment[0].length(); i++) {
		if (alignment[0][i] == alignment[1][i]) {
			count++;
		}
	}
	double pct = (double)count / (double)alignment[0].length();
	return pct;
}
void Graph::init(const double matrix[][4])
{
	graph = new double[l1*l2*3];
	backtrack = new tuple<int,int,int>[l1*l2*3];
	for (int i = 0; i < l1*l2*3; i++) {
		backtrack[i] = {-1,-1,-1};
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			score[i][j] = matrix[i][j];
		}
	}
	double highest = -1;
	for (int i = 0; i < l1; i++) {
		for (int j = 0; j < l2; j++) {
			fill_graph_at(i,j);
			int h = hash(i, j, 1);
			if (graph[h] > highest) {
				highest = graph[h];
				local_end = {i, j, 1};
			}
		}
	}
}

double Graph::get_score(char c1, char c2) const
{
	int t1 = transform(c1);
	int t2 = transform(c2);
	return score[t1][t2];
}

void Graph::fill_graph_at(int i, int j)
{
	if (i != 0) { // lower points
		auto p1 = graph[hash(i-1,j,0)] - epsilon;
		auto p2 = graph[hash(i-1,j,1)] - sigma;
		if (p1 > p2) {
			graph[hash(i,j,0)] = p1;
			backtrack[hash(i,j,0)] = {i-1, j, 0};
		} else {
			graph[hash(i,j,0)] = p2;
			backtrack[hash(i,j,0)] = {i-1, j, 1};
		}
	}
	if (j != 0) { // higher points
		//cout << "[" << i << "][" << j << "][2]: ";
		auto p1 = graph[hash(i,j-1,2)] - epsilon;
		auto p2 = graph[hash(i,j-1,1)] - sigma;
		if (p1 > p2) {
			//cout << "(" << i << ", " << j-1 << ", 2)" << endl;
			graph[hash(i,j,2)] = p1;
			backtrack[hash(i,j,2)] = {i, j-1, 2};
		} else {
			//cout << "(" << i << ", " << j-1 << ", 1)" << endl;
			graph[hash(i,j,2)] = p2;
			backtrack[hash(i,j,2)] = {i, j-1, 1};
		}
	}
	if (i == 0 && j != 0) { // lower points that depend on higher level
		graph[hash(i,j,0)] = graph[hash(i,j,2)];
		backtrack[hash(i,j,0)] = {i, j, 2};
	}
	if (j == 0 && i != 0) { // higher points that depend on lower level
		graph[hash(i,j,2)] = graph[hash(i,j,0)];
		backtrack[hash(i,j,2)] = {i, j, 0};
	}
	if (i != 0 || j != 0) { // middle level points not at origin

		vector<double> points = {graph[hash(i,j,0)],
					 graph[hash(i,j,2)]};
		if (i != 0 && j != 0) {
			points.push_back(graph[hash(i-1,j-1,1)] + get_score(s1[i-1], s2[j-1]));
		}
		if (local) {
			points.push_back(0);
		}
		double p_max = points[0];
		for (auto p : points) {
			if (p > p_max) {
				p_max = p;
			}
		}
		graph[hash(i,j,1)] = p_max;
//		cout << "[" << i << "][" << j << "][1]: ";
		if (points[0] == p_max) { // lower
			backtrack[hash(i,j,1)] = {i, j, 0};
//			cout << "(" << i << ", " << j << ", 0)" << endl;
		} else if (points[1] == p_max) {
			backtrack[hash(i,j,1)] = {i, j, 2};
//			cout << "(" << i << ", " << j << ", 2)" << endl;
		} else if (i != 0 && j != 0 && points[2] == p_max) {
			backtrack[hash(i,j,1)] = {i-1, j-1, 1};
//			cout << "(" << i-1 << ", " << j-1 << ", 1)" << endl;
		} else if (local) {
			backtrack[hash(i,j,1)] = {i, j, -2};
		}
	} else { // i == j == 0
		graph[hash(0,0,0)] = -100000;
		graph[hash(0,0,1)] = 0;
		graph[hash(0,0,2)] = -100000;
	}
}

vector<string> Graph::get_alignment(std::pair<int,int> &begin, std::pair<int,int> &end) const
{
	vector<string> alignment;
	vector<tuple<int,int,int>> path;
	tuple<int,int,int> pt = {l1-1, l2-1, 1};
	if (local) {
		pt = local_end;
	}
	int si1 = 0, si2 = 0;
	string str1, str2;
        while (get<2>(pt) != -1) {
		path.push_back(pt);
		if (local && get<2>(pt) == -2) {
			break;
		}
		si1 = get<0>(pt);
		si2 = get<1>(pt);
		auto next = backtrack[hash(get<0>(pt), get<1>(pt), get<2>(pt))];
		pt = next;
	}
	reverse(path.begin(), path.end());
	begin = std::make_pair(get<0>(path[0]),
			       get<1>(path[0]));
	end = std::make_pair(get<0>(path[path.size()-1]),
			     get<1>(path[path.size()-1]));
	for (int i = 1; i < path.size(); i++) {
		int x = get<0>(path[i]) - get<0>(path[i-1]);
		int y = get<1>(path[i]) - get<1>(path[i-1]);
//		int z = get<2>(path[i]) - get<2>(path[i-1]);
		if (x > 0) {
			str1 += s1[si1++];
			if (y == 0) {
				str2 += '-';
			}
		}
		if (y > 0) {
			str2 += s2[si2++];
			if (x == 0) {
				str1 += '-';
			}
		}
	}
	alignment.push_back(str1);
	alignment.push_back(str2);
	return alignment;
}
#endif
