/* -*- C++ -*-
 *
 * fastareader.h
 *
 * Author: Benjamin T James
 */

#ifndef FASTAREADER_H
#define FASTAREADER_H
#include "chrom.h"
#include <vector>
#include <functional>
class fastareader {
public:
	fastareader(std::vector<std::string> files_) : files(files_) {};
	void foreach_fasta(std::function<void(std::string, const std::string&)> func) const;
private:
	std::vector<std::string> files;
};
#endif
