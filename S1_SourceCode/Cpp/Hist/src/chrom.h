/* -*- C++ -*-
 *
 * chrom.h
 *
 * Author: Benjamin T James
 */

#ifndef CHROM_H
#define CHROM_H

#include <string>
class chrom {
public:
	virtual ~chrom() {};
	virtual size_t length() const = 0;
	virtual char at(size_t) const = 0;
	virtual std::string get_header() const = 0;
	virtual std::string substr(size_t, size_t) const = 0;
	virtual std::string get_base() const = 0;
};
#endif
