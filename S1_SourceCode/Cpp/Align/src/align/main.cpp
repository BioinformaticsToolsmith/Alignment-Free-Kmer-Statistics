#include <cmath>
#include <map>
#include <mutex>
#include <thread>
#include "DNA.hpp"
#include "Graph.hpp"
#ifdef HIST
#include "../utility/Util.h"
#include "../nonltr/KmerHashTable.h"
#include "../nonltr/ChromosomeOneDigit.h"
#include "../nonltr/ChromListMaker.h"
#endif
using namespace std;

#ifdef OPENMP
#include <omp.h>
#endif

// const double global_matrix[4][4] = {{1, -3, -3, -3},
// 				    {-3, 1, -3, -3},
// 				    {-3, -3, 1, -3},
// 				    {-3, -3, -3, 1}};
// const double global_sigma = 5;
// const double global_epsilon = 2;
const double global_matrix[4][4] = {{1, -1, -1, -1},
				    {-1, 1, -1, -1},
				    {-1, -1, 1, -1},
				    {-1, -1, -1, 1}};
const double global_sigma = 1;
const double global_epsilon = 1;

bool isInside(const std::string & str, char c)
{
    return str.find(c) != std::string::npos;
}

void infer(int x, int y, int n1[4], int n2[4], int mat[4][4])
{
	mat[x][y]++;
	n1[x]++;
	n2[y]++;
}
int comp(char a) {
	int x = transform(a);
	return ((x + 2) % 4);
}
void observe(char x, char y, int n1[4], int n2[4], int mat[4][4])
{
	if (x == '-' && y != '-') {
		//n2[transform(y)]++;
	} else if (y == '-' && x != '-') {
		//n1[transform(x)]++;
	} else if (x != '-' && y != '-') {
		infer(transform(x), transform(y), n1, n2, mat);
		// infer(comp(x), comp(y), n1, n2, mat);
		infer(transform(y), transform(x), n1, n2, mat);
		// infer(comp(y), comp(x), n1, n2, mat);
	}
}

void get_mat_data(const vector<const DNA*> &vec, double &sigma, double &epsilon, double (&s)[4][4])
{
	int n1[4] = {4, 4, 4, 4}, n2[4] = {4, 4, 4, 4};
	int count = 0;
	int mat[4][4] = {{1, 1, 1, 1},
			 {1, 1, 1, 1},
			 {1, 1, 1, 1},
			 {1, 1, 1, 1}};
	for (int i = 0; i < vec.size(); i++) {
		for (int j = i + 1; j < vec.size(); j++) {
			Graph g(vec[i]->get_chrom(), vec[j]->get_chrom(), global_matrix, global_sigma, global_epsilon, true);
			std::pair<int,int> begin, end;
			auto alignment = g.get_alignment(begin, end);
			// if (isInside(alignment[0], '-') || isInside(alignment[1], '-')) {
			// 	continue;
			// }
			double identity = g.get_identity(alignment);
			//cout << alignment[0] << " " << alignment[1] << " " << identity << endl;
			// if (identity > 0.7) {
			// 	continue;
			// }
			// if (alignment[0].length() > 0.7 * vec[i]->get_chrom().length()) {
			// 	cout << "identity" << endl;
			// 	continue;
			// }
			count++;
			if (count == 100) {
				i = vec.size();
				break;
			} else {
				int strl = alignment[0].length();
				for (int i = 0; i < strl; i++) {
					observe(alignment[0][i], alignment[1][i], n1, n2, mat);
				}
			}
		}
	}
	double q1[4]; // normalized n1
	double q2[4]; // normalized n2
	double p[4][4];              // normalized mat
	int npairs = n1[0] + n1[1] + n1[2] + n1[3];
	for (int x = 0; x < 4; x++) {
		q1[x] = (double)n1[x] / npairs;
		q2[x] = (double)n2[x] / npairs;
		for (int y = 0; y < 4; y++) {
			p[x][y] = (double)mat[x][y] / npairs;
		}
	}
        double largest = 0;
	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			s[x][y] = log(p[x][y] / (q1[x] * q2[y]));
			if ((x == 0 && y == 0) || s[x][y] > largest) {
				largest = s[x][y];
			}
		}
	}
	count = 0;
	double sum = 0;
	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			int i = 100.0 * s[x][y] / largest;
			double d = (double)i / 100.0;
			if (x != y) {
				sum += d;
				count++;
			}
			//cout << i << "  ";
			printf("%0.2f  ", d);
//			cout << d << "  ";
		}
		cout << endl;
	}
	double avg = sum / count;
	sigma = fabs(avg / 2.0);
	epsilon = sigma * 0.75;
	cout << "sigma: " << sigma << endl;
	cout << "epsilon: " << epsilon << endl;

}

#ifdef HIST
void get_histograms(int k, vector<string> fileList, string out)
{
	ofstream fout;
	fout.open(out);
	for (auto file : fileList) {
		ChromListMaker *maker = new ChromListMaker(file);
		auto chromList = maker->makeChromOneDigitList();
		for (auto chm : *chromList) {
			auto chrom = dynamic_cast<ChromosomeOneDigit*>(chm);
			KmerHashTable<unsigned long, double> table(k, 1);
			const vector<vector<int>*> *segment = chrom->getSegment();
			vector<int> kmers;
			const char *seg_bases = chrom->getBase()->c_str();
			int len = chrom->size();
			int i = 0;
			for (auto s : *segment) {
				int start = s->at(0);
				int end = s->at(1);
				table.wholesaleIncrement(seg_bases, start, end - k + 1);
				// if ((3 * i) % len == 0) {
				// 	kmers.push_back(table.valueOf(seg_bases, i));
				// }
				// i++;
			}
			vector<const char*> keys;
			table.getKeys(keys);
			fout << chrom->getHeader();
			for (const char *str : keys) {
				fout << "," << table.valueOf(str);
			}
			fout << endl;
		}
	}
	fout.close();
}
#endif

//mutex hashmap_lock;

// void do_alignment(int i, const vector<const DNA*> vec, const double mat[4][4], double sig, double eps)
// {
// 	for (int j = i; j < vec.size(); j++) {
// 		Graph g(vec[i]->get_chrom(), vec[j]->get_chrom(), mat, sig, eps, false);
// 		auto score = g.get_total_score();
// 		auto identity = g.get_identity(g.get_alignment());
// 		auto pr = make_pair(make_pair(i, j), make_pair(score, identity));
// 		//unlock
// //		lock_guard<mutex> guard(hashmap_lock);
// 		hashmap.insert(pr);
// 		//lock
// 	}
// }

void usage(string progname)
{
	cerr << "Usage: " << progname << " file1.fa file2.fa dir1/*.fa -k 4 -o output.hist" << endl;
	cerr << "    or " << progname << " file1.fa file2.fa dir1/*.fa -d , -o output.align" << endl;
	cerr << "    or " << progname << " file1.fa file2.fa dir1/*.fa -l -d , -o output.localalign" << endl;
	cerr << endl;
	cerr << "Output column format:" << endl;
	cerr << "First column:      >header 1" << endl;
	cerr << "Second column:     >header 2" << endl;
	cerr << "Third column:      Internal alignment score [0-1]" << endl;
	cerr << "Fourth column:     Alignment Identity score [0-1]" << endl;
	cerr << "Fifth column:      Coverage score [0-1]" << endl;
	cerr << "Sixth column:      seq1 index (for sorting the output)" << endl;
	cerr << "Seventh column:    seq2 index (for sorting the output)" << endl;
	cerr << endl;
	cerr << "For making this output compatible with all 3 experiments, only 4 columns are used:" << endl;
	cerr << "First column, second column, third column, and an alignment score for the 4th column." << endl;
	cerr << "Also, the data has to be sorted first by column 6, then by column 7 to match histogram files." << endl;
	cerr << endl;


}
void get_args(vector<string> args, vector<string>& files, string& outputfile, char& delim, int& k, bool& is_local)
{
	for (int i = 1; i < args.size(); i++) {
		if (args[i][0] == '-' && i + 1 < args.size() && args[i].length() == 2) {
			switch (args[i][1]) {
			case 'k':
			case 'K':
				k = stod(args[++i]);
				if (k <= 0) {
					usage(args[0]);
					exit(1);
				}
				break;
			case 'd':
			case 'D':
				delim = args[++i][0];
				break;
			case 'l':
			case 'L':
				is_local = true;
				break;
			case 'o':
			case 'O':
				outputfile = args[++i];
				break;
			case 'h':
				usage(args[0]);
				exit(0);
			default:
				cerr << "bad option \"" << args[i] << "\"" << endl;
				usage(args[0]);
				exit(1);
			}
		} else {
			files.push_back(args[i]);
		}
	}
	if (outputfile == "" || files.empty()) {
		usage(args[0]);
		exit(1);
	}
}
// double align(string sa, string sb)
// {
// 	int la = sa.length();
// 	int lb = sb.length();
// 	utility::LCSLen lcs(sa.c_str(), 0, la - 1,
// 			    sb.c_str(), 0, lb - 1);
// 	int sublen = lcs.getLenCS();
// 	int denom = std::max(la, lb);
// 	double sim = 1.0 * sublen / denom;
// 	return sim;
// }
int main(int argc, char *argv[])
{
	srand(time(NULL));
	vector<string> args;
	for (int i = 0; i < argc; i++) { args.push_back(argv[i]); }

	int k = -1;
	char delim = ';';
	vector<string> files;
	string output_file;
	bool is_local = false;
	get_args(args, files, output_file, delim, k, is_local);
	if (k != -1) {
#ifdef HIST
		get_histograms(k, files, output_file);
#endif
		return 0;
	}




	auto vec = get_sequences(files);
	double mat[4][4];
	double sigma, epsilon;
	sigma = global_sigma;
	epsilon = global_epsilon;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = global_matrix[i][j];
		}
	}

//	get_mat_data(vec, sigma, epsilon, mat);

	// int count = 0;
	// for (int i = 0; i < vec.size(); i++) {
	// 	for (int j = i; j < vec.size(); j++) {
	// 		count++;
	// 	}
	// }
	// int alignments = 0;
	// for (int i = 0; i < vec.size(); i += num_threads) {
	// 	vector<thread*> threads;
	// 	for (int t = 0; t < num_threads && i+t < vec.size(); t++) {
	// 		threads.push_back(new thread(do_alignment, i+t, vec, mat, sigma, epsilon));
	// 		alignments += vec.size() - (i + t);
	// 	}
	// 	for (auto t : threads) {
	// 		t->join();
	// 		delete t;
	// 	}
	// }
	// #ifdef OPENMP
        // #pragma omp parallel for shared(hashmap)
	// #endif
	int total_size = vec.size() * (vec.size() + 1) / 2;
	unsigned long count = 0;
	double last_pct = 0;
	ofstream f;
	f.open(output_file);
	#ifdef OPENMP
#pragma omp parallel for schedule(dynamic)
	#endif
	for (int i = 0; i < vec.size(); i++) {
		for (int j = i; j < vec.size(); j++) {
			double identity;
			int ilen = vec[i]->get_chrom().length();
			int jlen = vec[j]->get_chrom().length();
			Graph g(vec[i]->get_chrom(), vec[j]->get_chrom(), mat, sigma, epsilon, is_local);
			auto score = g.get_total_score();
			std::pair<int,int> begin, end;
			auto algnment = g.get_alignment(begin, end);
			double algnlen = 0.5 * (algnment[0].length() + algnment[1].length());
			identity = g.get_identity(algnment);

			auto pr = make_pair(make_pair(i, j), make_pair(score, identity));
			double coverage = algnlen / (std::max(begin.first, begin.second)
						     + algnlen
						     + std::max(ilen - end.first,
								jlen - end.second));
#pragma omp critical
			{

				count++;
				//hashmap.insert(pr);
				f << vec[i]->get_header() << delim;
				f << vec[j]->get_header() << delim;
				f << pr.second.first << delim << pr.second.second << delim;
				// f << coverage << delim;
				// f << ilen << delim << jlen << delim;
				// f << begin.first << delim << end.first << delim;
				// f << begin.second << delim << end.second << delim;
				f << i << delim << j << endl;
				f.flush();
				double pct = 100.0 * count / total_size;
				if (pct - last_pct > 0.2) {
					int amt = round(pct * 100.0);
					cout << amt / 100.0 << "%" << endl;
					last_pct = pct;
				}
			}
		}
	}
/*

	for (int i = 0; i < vec.size(); i++) {
		for (int j = i + 1; j < vec.size(); j++) {
			Graph g(vec[i]->get_chrom(), vec[j]->get_chrom(), sigma, epsilon, false);
			f << vec[i]->get_header() << "," << vec[j]->get_header() << ",";

			auto score = g.get_total_score();
			f << score << ",";
			auto identity = g.get_identity(g.get_alignment());
			hashmap.insert(make_pair(make_pair(i, j), make_pair(score, identity)));
			f << identity << endl;
		}
	}
*/
	// ofstream f;
	// f.open(output_file);
	// for (int i = 0; i < vec.size(); i++) {
	// 	for (int j = i; j < vec.size(); j++) {
	// 		f << vec[i]->get_header() << delim;
	// 		f << vec[j]->get_header() << delim;
	// 		auto pr = hashmap[make_pair(i, j)];
	// 		f << pr.first << delim << pr.second << endl;
	// 	}
	// }
	f.close();
	return 0;
}
