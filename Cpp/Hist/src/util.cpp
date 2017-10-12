/* -*- C++ -*-
 *
 * util.cpp
 *
 * Author: Benjamin T James
 */

#include "util.h"
#include <cmath>
#include <iostream>
int util::gcd(int a, int b)
{
	if (b <= 0) {
		return a;
	}
	return gcd(b, a % b);
}

char util::decode(uint8_t u)
{
	char c;
	switch (u) {
	case 0:
		c = 'A';
		break;
	case 1:
		c = 'C';
		break;
	case 2:
		c = 'G';
		break;
	case 3:
		c = 'T';
		break;
	default:
		throw "bad";
	}
	return c;
}
uint8_t util::encode(char c)
{
	uint8_t ret = 0;
	int r;
	switch (c) {
	case 'a':
	case 'A':
		ret = 0;
		break;
	case 'c':
	case 'C':
		ret = 1;
		break;
	case 'g':
	case 'G':
		ret = 2;
		break;
	case 't':
	case 'T':
	case 'u':
	case 'U':
		ret = 3;
		break;
	case 'r':
	case 'R': // a or g
		if (rand() % 2) {
			ret = encode('A');
		} else {
			ret = encode('G');
		}
		break;
	case 'y':
	case 'Y': // c, t, or u
		r = rand() % 2;
		if (r == 0) {
			ret = encode('C');
		} else {
			ret = encode('T');
		}
		break;
	case 'k':
	case 'K': // g, t, or u
		r = rand() % 2;
		if (r == 0) {
			ret = encode('G');
		} else {
			ret = encode('T');
		}
		break;
	case 'm':
	case 'M':
		r = rand() % 2;
		if (r == 0) {
			ret = encode('A');
		} else {
			ret = encode('C');
		}
		break;
	case 's':
	case 'S':
		r = rand() % 2;
		if (r == 0) {
			ret = encode('C');
		} else {
			ret = encode('G');
		}
		break;
	case 'w':
	case 'W':
		r = rand() % 2;
		if (r == 0) {
			ret = encode('A');
		} else {
			ret = encode('T');
		}
		break;
	case 'b':
	case 'B':
		r = rand() % 3;
		if (r == 0) {
			ret = encode('C');
		} else if (r == 1) {
			ret = encode('G');
		} else {
			ret = encode('T');
		}
		break;
	case 'd':
	case 'D':
		r = rand() % 3;
		if (r == 0) {
			ret = encode('A');
		} else if (r == 1) {
			ret = encode('G');
		} else {
			ret = encode('T');
		}
		break;
	case 'h':
	case 'H':
		r = rand() % 3;
		if (r == 0) {
			ret = encode('A');
		} else if (r == 1) {
			ret = encode('C');
		} else {
			ret = encode('T');
		}
		break;
	case 'v':
	case 'V':
		r = rand() % 3;
		if (r == 0) {
			ret = encode('A');
		} else if (r == 1) {
			ret = encode('C');
		} else {
			ret = encode('G');
		}
		break;
	case 'n':
	case 'N':
		ret = rand() % 4;
		break;
	default:
		throw "bad nucl";
	}
	return ret;
}

int util::gcd_vec(std::vector<int> v)
{
	int ret = v[0];
	for (size_t i = 1; i < v.size(); i++) {
		if (v[i] == 0) {
			continue;
		}
		ret = gcd(ret, v[i]);
	}
	return ret;
}

inline int sign(double x) {
	return (x > 0) - (x < 0);
}

void util::scale(double (&mat)[4][4], double &sigma, double& epsilon)
{
	double scale_factor = 100000;
	std::vector<int> signs, scaled;
	signs.push_back(sign(sigma));
	scaled.push_back(round(scale_factor * fabs(sigma)));
	signs.push_back(sign(epsilon));
	scaled.push_back(round(scale_factor * fabs(epsilon)));
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			signs.push_back(sign(mat[i][j]));
			scaled.push_back(round(scale_factor * fabs(mat[i][j])));
		}
	}
	double common_div = gcd_vec(scaled);
	sigma = signs[0] * scaled[0] / common_div;
	epsilon = signs[1] * scaled[1] / common_div;
	int count = 2;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = signs[count] * scaled[count] / common_div;
			count++;
		}
	}
}
