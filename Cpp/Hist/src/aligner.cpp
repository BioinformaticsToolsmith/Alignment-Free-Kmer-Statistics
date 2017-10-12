/* -*- C++ -*-
 *
 * aligner.cpp
 *
 * Author: Benjamin T James
 */

#include "aligner.h"
#include "needleman_wunsch.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <sstream>
void aligner::doall(double (&dmatrix)[4][4], double dsigma, double depsilon)
{
	std::ofstream out(output);
	util::scale(dmatrix, dsigma, depsilon);
	uintmax_t num_completed = 0;
	uintmax_t num_to_do = vec.size() * (vec.size() + 1) / 2;
	double ratio = 0;
	double lastratio = 0;
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < vec.size(); i++) {
		std::string si = vec[i].get_base();
		std::string hi = vec[i].get_header();
		for (size_t j = i; j < vec.size(); j++) {
			std::string sj = vec[j].get_base();
			std::string hj = vec[j].get_header();
			needleman_wunsch nw(si, sj, dmatrix, dsigma, depsilon);
			std::pair<std::string, std::string> algn;
			if (i == j) {
				algn = make_pair(si, sj);
			} else {
				algn = nw.align();
			}
			double sim = nw.identity(algn);
#pragma omp critical
			{
				ratio = 100.0 * ++num_completed / num_to_do;
				if (ratio - lastratio >= 0.1) {
					std::cout << "\rProgress: " << ratio << "%   ";
					std::cout.flush();
					lastratio = ratio;
				}
				if (i != j) {
					output_sim(out, hj, hi, sim);
				}
				output_sim(out, hi, hj, sim);
			}
		}
	}
	std::cout << std::endl;
}


void aligner::output_sim(std::ofstream &out, std::string h1, std::string h2, double sim) const
{
	out << h1 << delim << h2 << delim << sim << std::endl;
}
