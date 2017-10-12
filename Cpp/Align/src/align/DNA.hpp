#ifndef DNA_HPP
#define DNA_HPP

#include <iostream>
#include <vector>
#include <dirent.h>
#include <fstream>
#include <regex>
#include "../nonltr/ChromosomeOneDigit.h"
#include "../nonltr/ChromListMaker.h"
using namespace std;

class DNA {
public:
	DNA(const string hdr, const string chr) : header(hdr), chrom(chr) {};
	string get_header() const { return header; };
	string get_chrom() const { return chrom; };
private:
	const string header, chrom;
};


vector<const DNA*> get_sequences(vector<string> filenames)
{
	vector<const DNA*> sequences;
	for (auto filename : filenames) {
		// fstream in(filename);
		// string line, chrom = "";
		// string header;
		// while (getline(in, line)) {
		// 	if (line[0] == '>') {
		// 		if (chrom.length() > 0) {
		// 			sequences.push_back(new DNA(header, chrom));
		// 			chrom = "";
		// 		}
		// 		header = line;
		// 	} else {
		// 		chrom += line;
		// 	}
		// }
		// if (chrom.length() > 0) {
		// 	sequences.push_back(new DNA(header, chrom));
		// }
		ChromListMaker *maker = new ChromListMaker(filename);
		auto chromList = maker->makeChromOneDigitList();
		for (auto chm : *chromList) {
			auto chrom = dynamic_cast<ChromosomeOneDigit*>(chm);
			sequences.push_back(new DNA(chrom->getHeader(), *chrom->getBase()));
		}
	}
	cout << "Got sequences" << endl;
	return sequences;
}
#endif
