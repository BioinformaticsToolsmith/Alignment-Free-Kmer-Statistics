/* -*- C++ -*-
 *
 * aligner.h
 *
 * Author: Benjamin T James
 */
#ifndef ALIGNER_H
#define ALIGNER_H

#include "bitchrom.h"
#include <vector>
#include <fstream>
class aligner {
public:
	aligner(const std::vector<bitchrom<int> > &vec_,
		std::string output_, std::string delim_)
		: vec(vec_), output(output_), delim(delim_) {};
	void doall(double (&mat)[4][4], double sigma, double epsilon);
private:
	void output_sim(std::ofstream &out, std::string h1, std::string h2, double sim) const;
	const std::vector<bitchrom<int> > &vec;
	std::string output, delim;
};
#endif
