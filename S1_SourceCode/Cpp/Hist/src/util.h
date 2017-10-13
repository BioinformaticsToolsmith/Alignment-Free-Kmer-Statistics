/* -*- C++ -*-
 *
 * util.h
 *
 * Author: Benjamin T James
 */
#ifndef UTIL_H
#define UTIL_H
#include <vector>
#include <stdint.h>
class util {
public:
	static void scale(double (&mat)[4][4], double &sigma, double& epsilon);
	static int gcd(int a, int b);
	static int gcd_vec(std::vector<int> v);
	static uint8_t encode(char c);
	static char decode(uint8_t);
};

#endif
