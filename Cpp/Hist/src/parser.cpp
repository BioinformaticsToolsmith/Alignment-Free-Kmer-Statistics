/* -*- C++ -*-
 *
 * parser.cpp
 *
 * Author: Benjamin T James
 */

#include "parser.h"
#include <iostream>
#include <fstream>

int input_int(std::string arg, int& r)
{
	try {
		int i = std::stoi(arg);
		if (i > 0) {
			r = i;
		} else {
			return -1;
		}
	} catch (const std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
		return -1;
	}
	return 0;
}

int read_matrix(std::string file, double (&mat)[4][4] )
{
	std::ifstream in(file);
	if (!in) {
		std::cerr << "Cannot open file \"" << file << "\"" << std::endl;
		return -1;
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (!in) {
				return -1;
			}
			in >> mat[i][j];
		}
	}
	return 0;
}
int input_double(std::string arg, double& r)
{
	try {
		r = std::stod(arg);
	} catch (const std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
		return -1;
	}
	return 0;
}

void parser::defaults(std::string argv0)
{
	prog_name = argv0;
	k = -1;
	delim = ",";
	output = "";
	sigma = -3;
	epsilon = -1;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				matrix[i][j] = 1;
			} else {
				matrix[i][j] = -1;
			}
		}
	}
	srand(time(NULL));

}
void parser::parse(int argc, char **argv)
{
	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			std::string arg = argv[i];
			if (arg.length() != 2 || i + 1 == argc) {
				usage(1);
			}
			switch (argv[i][1]) {
			case 'h':
				usage();
				break;
			case 's':
				if (input_double(argv[++i], sigma) == -1) {
					usage(1);
				}
				break;
			case 'e':
				if (input_double(argv[++i], epsilon) == -1) {
					usage(1);
				}
				break;
			case 'm':
				if (read_matrix(argv[++i], matrix) == -1) {
					usage(1);
				}
				break;
			case 'k':
				if (input_int(argv[++i], k) == -1) {
					usage(1);
				}
				break;
			case 'd':
				delim = argv[++i];
				break;
			case 'o':
				output = argv[++i];
				break;
			}
		} else {
			files.push_back(argv[i]);
		}
	}
	if (files.empty() || output.length() == 0) {
		usage(1);
	}
	mode = (k == -1) ? "align" : "hist";

}
parser::parser(int argc, char **argv)
{
	defaults(*argv);
	parse(argc, argv);
}

void parser::usage(int exitval) const
{
	std::cout << "Usage for histograms:" << std::endl;
	std::cout << "\t" << prog_name;
	std::cout << " -k 5 *.fasta -o output [-d delim]" << std::endl;
	std::cout << "Usage for all vs all alignment:" << std::endl;
	std::cout << "\t" << prog_name;
	std::cout << " *.fasta -o output [-d delim] [-s sigma] [-e epsilon] [-m matrixfile]" << std::endl;
	std::cout << std::endl;
	std::cout << "Matrix file is 4 lines with 4 numbers on each line" << std::endl;
	std::cout << "separated by spaces, by default 1 for match -1 for mismatch" << std::endl;
	std::cout << std::endl;
	std::cout << "Sigma: affine gap opening penalty, default -3" << std::endl;
	std::cout << "Epsilon: affine gap extension cost, default -1" << std::endl;
	exit(exitval);
}
