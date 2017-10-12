/* -*- C++ -*-
 *
 * kmercounter.cpp
 *
 * Author: Benjamin T James
 */

#include "kmercounter.h"
#include <vector>
#include <cmath>
#include "util.h"
kmercounter::kmercounter(int k_, std::string output_file, std::string delim_)
	: k(k_), delim(delim_)
{
	out.open(output_file);
}

size_t hash(const std::string& str, ssize_t idx, int k)
{
	size_t sum = 0;
	for (ssize_t i = idx; i < idx + k; i++) {
		sum = 4 * sum + util::encode(str[i]);
	}
	return sum;
}
std::vector<int> get_histogram(const std::string& str, int k)
{
	int four_k = ceil(pow(4, k));
	std::vector<int> vec(four_k);
	std::fill_n(vec.begin(), vec.size(), 1);
	const ssize_t len = str.length();
	for (ssize_t i = 0; i + k <= len; i++) {
		size_t idx = hash(str, i, k);
		vec[idx]++;
	}
	return vec;
}
void kmercounter::count(std::string header, const std::string& body)
{
	auto hist = get_histogram(body, k);
	out << header;
	for (auto kmer : hist) {
		out << delim << kmer;
	}
	out << std::endl;
}
