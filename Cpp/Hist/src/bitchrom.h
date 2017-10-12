/* -*- C++ -*-
 *
 * bitchrom.h
 *
 * Author: Benjamin T James
 */
#ifndef BITCHROM_H
#define BITCHROM_H

#include "chrom.h"
#include <vector>

template<class T>
class bitchrom : public chrom {
public:
	bitchrom(std::string header_, const std::string & base);
	size_t length() const { return len; };
	char at(size_t index) const { return substr(index, 1)[0]; };
	std::string get_header() const { return header; };
	std::string substr(size_t begin, size_t slen) const;
	std::string get_base() const { return substr(0, len); };
private:
	uint8_t encode(char c) const;
	char decode(uint8_t) const;
	void fill(const std::string & base);
	std::vector<T> data;
	const std::string header;
	const size_t len;
	const size_t letters_per_datum;
};
#include "bitchrom.hpp"
#endif
