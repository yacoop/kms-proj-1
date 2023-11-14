#pragma once

#define _USE_MATH_DEFINES
#include "vec3.h"
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

struct Parameters {
	double m = 0;
	double e = 0;
	int n = 0;
	double R = 0;
	double f = 0;
	double a = 0;
	double T_0 = 0;
	double tau = 0;
	int S_o = 0;
	int S_d = 0;
	int S_out = 0;
	int S_xyz = 0;
	double k = 0;
	int N = 0;
	double L = 0;

	Parameters(const char* filename) {
		//initiate parameters from the file
		std::ifstream file(filename);

		std::string s;
		file >> n;
		std::getline(file, s);
		file >> m;
		std::getline(file, s);
		file >> e;
		std::getline(file, s);
		file >> R;
		std::getline(file, s);
		file >> f;
		std::getline(file, s);
		file >> L;
		std::getline(file, s);
		file >> a;
		std::getline(file, s);
		file >> T_0;
		std::getline(file, s);
		file >> tau;
		std::getline(file, s);
		file >> S_o;
		std::getline(file, s);
		file >> S_d;
		std::getline(file, s);
		file >> S_out;
		std::getline(file, s);
		file >> S_xyz;
		std::getline(file, s);

		N = n * n * n;
		L = 1.22 * a * ((double)n - 1) + 0.1;
		k = 8.31e-3;
	}
};

struct State {
	double t = 0;
	double E = 0;
	double V = 0;
	double H = 0;
	double T = 0;
	double P = 0;


	std::vector<vec3> r;
	std::vector<vec3> p;
	std::vector<vec3> F;
};