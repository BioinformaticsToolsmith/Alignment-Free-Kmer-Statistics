/* -*- C++ -*-
 *
 * main.cpp
 *
 * Author: Benjamin T James
 */

#include "fastareader.h"
#include "parser.h"
#include "kmercounter.h"
#include "aligner.h"
#include "bitchrom.h"
#include <iostream>

int main(int argc, char **argv)
{
	parser p(argc, argv);
	auto files = p.get_files();
	auto output_file = p.get_output();
	auto delim = p.get_delim();

	if (p.get_mode() == "align") {
		fastareader reader(files);
		std::vector<bitchrom<int> > vec;
		reader.foreach_fasta([&](std::string header, const std::string& body) {
				vec.emplace_back(header, body);
			});
		aligner aln(vec, output_file, delim);
		double matrix[4][4];
		p.get_matrix(matrix);
		aln.doall(matrix, p.get_sigma(), p.get_epsilon());
	} else {
		fastareader reader(files);
		int k = p.get_k();
		kmercounter cnt(k, output_file, delim);
		reader.foreach_fasta([&](std::string header, const std::string& body) {
				cnt.count(header, body);
			});
	}
	return 0;
}
