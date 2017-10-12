/* -*- C++ -*-
 *
 * kmercounter.h
 *
 * Author: Benjamin T James
 */

#ifndef KMERCOUNTER_H
#define KMERCOUNTER_H
#include <string>
#include "chrom.h"
#include <fstream>
class kmercounter {
public:
	kmercounter(int k_, std::string output_file_, std::string delim_);
	void count(std::string header, const std::string& body);
private:
	const int k;
	std::ofstream out;
	std::string delim;
};
#endif
