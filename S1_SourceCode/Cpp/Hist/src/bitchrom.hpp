/* -*- C++ -*-
 *
 * bitchrom.hpp
 *
 * Author: Benjamin T James
 */

#include "bitchrom.h"
#include "util.h"
#include <cmath>


template<class T>
bitchrom<T>::bitchrom(std::string header_, const std::string& base)
	: header(header_), len(base.length()), letters_per_datum(4 * sizeof(T))
	// each byte can store 4 letters (2 bit integers)
{
	data.reserve(1 + ceil((double)len / letters_per_datum));
	for (size_t i = 0; i < len; i++) {
		data.push_back(0); // there's probably a better way to do this
	}
	fill(base);
}

template<class T>
void bitchrom<T>::fill(const std::string &base)
{
	for (size_t i = 0; i < len; i++) {
		size_t offset = 2 * (i % letters_per_datum);
		size_t block = i / letters_per_datum;
		data[block] |= (util::encode(base[i]) & 0x3) << offset;
	}
}

template<class T>
std::string bitchrom<T>::substr(size_t begin, size_t slen) const
{
	if (begin + slen > len) {
		throw "bad length";
	}
	size_t last_idx = begin + slen - 1;
	std::string s = "";
	for (size_t i = begin; i <= last_idx; i++) {
		size_t offset = 2 * (i % letters_per_datum);
		size_t block = i / letters_per_datum;
		uint8_t encoded = (data[block] >> offset) & 0x3;
		s += util::decode(encoded);
	}
	return s;
}
