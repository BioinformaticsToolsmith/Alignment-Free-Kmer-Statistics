/* -*- C++ -*-
 *
 * fastareader.cpp
 *
 * Author: Benjamin T James
 */

#include "fastareader.h"
#include <fstream>
void foreach_file(std::string filename, std::function<void(std::string header, const std::string& body)> func)
{
	std::ifstream in(filename);
	std::string header = "", body = "";
	for (std::string line; getline(in, line);) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '>') {
			if (!header.empty()) {
				func(header, body);
			}
			header = line;
			body.clear();
		} else {
			body += line;
		}
	}
	if (!header.empty()) {
		func(header, body);
	}
}
void fastareader::foreach_fasta(std::function<void(std::string header, const std::string& body)> func) const
{
	for (std::string fname : files) {
		foreach_file(fname, func);
	}
}
