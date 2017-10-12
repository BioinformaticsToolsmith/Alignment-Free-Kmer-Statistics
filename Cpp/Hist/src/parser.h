/* -*- C++ -*-
 *
 * parser.h
 *
 * Author: Benjamin T James
 */

#ifndef PARSER_H
#define PARSER_H
#include <vector>
#include <string>

class parser {
public:
	parser(int argc, char **argv);
	void usage(int exitval=0) const;
	std::vector<std::string> get_files() const { return files; }
	std::string get_mode() const { return mode; }
	int get_k() const { return k; }
	std::string get_delim() const { return delim; }
	std::string get_output() const { return output; }
	void get_matrix(double (&mat)[4][4]) const {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				mat[i][j] = matrix[i][j];
			}
		}
	}
	double get_sigma() const { return sigma; }
	double get_epsilon() const { return epsilon; }
private:
	void defaults(std::string argv0);
	void parse(int argc, char **argv);
	std::vector<std::string> files;
	std::string prog_name, mode, delim, output;
	int k;
	double sigma, epsilon;
	double matrix[4][4];
};
#endif
