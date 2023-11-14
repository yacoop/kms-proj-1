#include "argon.h"

class Simulation {
public:
	Simulation(const char* parametersFile) : params(parametersFile) {
		InitState();
	}

	void run() {
		//seed random
		srand((unsigned)time(0));
		//open files
		std::ofstream xyzFile("output.xyz");
		std::ofstream stateFile("output.csv");

		std::vector<vec3> p_temp;

		//run simulation
		for (int i = 0; i < params.S_d; i++) {

			p_temp = state.p;
			for (int i = 0; i < params.N; i++) {
				p_temp[i] += 0.5 * state.F[i] * params.tau;
			}
			for (int i = 0; i < params.N; i++) {
				state.r[i] += p_temp[i] * params.tau / params.m;
			}

			calc_FVP();
			
			for (int i = 0; i < params.N; i++) {
				state.p[i] = p_temp[i] + 0.5 * state.F[i] * params.tau;
			}
			
			calc_EnT();
	
			state.t += params.tau;

			if (i % params.S_xyz == 0)
				print_xyz(xyzFile);

			if (i % params.S_out == 0)
				print_state(stateFile);
		}



		//close files
		xyzFile.close();
		stateFile.close();
	}

private:
	void InitState() {
		state.r.resize(params.N);
		state.p.resize(params.N);
		state.F.resize(params.N);

		double a = params.a;
		int n = params.n;

		//initiate r
		vec3 b_0{ a, 0, 0 };
		vec3 b_1{ a / 2, a * std::sqrt(3.0) / 2, 0 };
		vec3 b_2{ a / 2, a * std::sqrt(3.0) / 6, a * std::sqrt(2.0 / 3) };
		vec3 r_vec;
		for (double i_0 = 0; i_0 < n; i_0++) {
			for (double i_1 = 0; i_1 < n; i_1++) {
				for (double i_2 = 0; i_2 < n; i_2++) {
					int i = (int)(i_0 + i_1 * n + i_2 * n * n);
					r_vec = (i_0 - (double)(n - 1) / 2) * b_0 + (i_1 - (double)(n - 1) / 2) * b_1 + (i_2 - (double)(n - 1) / 2) * b_2;
					state.r[i] = r_vec;
				}
			}
		}

		//initiate p
		vec3 sum_p = vec3(0, 0, 0);
		for (int i = 0; i < params.N; i++) {
			state.p[i] = { rand_momentum(), rand_momentum(), rand_momentum() };
			sum_p += state.p[i];
		}

		for (int i = 0; i < params.N; i++) {
			state.p[i] -= sum_p / params.N * vec3(1, 1, 1);
		}

		
		calc_FVP();
		calc_EnT();
	}

	double rand_momentum()
	{
		double r = rand() % 2;
		double s;
		if (r == 0)
			s = -1;
		else
			s = 1;

		r = (double) rand() / RAND_MAX;
		double E = - 0.5 * params.k * params.T_0 * std::log(r);
		return s * std::sqrt(2 * params.m * E);
	}


	void calc_FVP() {
		state.F.clear();
		state.F.resize(params.N);

		state.V = 0;
		state.P = 0;

		double r_i, r_ij2, r_part6, r_part12;
		vec3 F_ij;

		for (int i = 0; i < params.N; i++) {
			r_i = state.r[i].length();

			if (r_i > params.L) {
				vec3 F_s = params.f * (params.L - r_i) * (state.r[i] / r_i);
				state.P += F_s.length();
				state.F[i] += F_s;
				state.V += 0.5 * params.f * (r_i - params.L) * (r_i - params.L);
			}
			for (int j = 0; j < i; j++) {
				r_ij2 = (state.r[i] - state.r[j]).length_squared();
				r_part6 = params.R*params.R/r_ij2;
				r_part6 = r_part6 * r_part6 * r_part6;
				r_part12 = r_part6 * r_part6;
				F_ij = 12 * params.e * (r_part12 - r_part6) * (state.r[i] - state.r[j]) / r_ij2;
				state.F[i] += F_ij;
				state.F[j] -= F_ij;
				state.V += params.e * (r_part12 - 2 * r_part6);
			}
		}

		state.H = state.V + state.E;
		state.P = state.P * 1.0/(4 * M_PI * params.L);
	}

	void calc_EnT()
	{
		state.E = 0;
		for(vec3 &p_i: state.p)
			state.E += p_i.length_squared() / (2 * params.m);
		state.T = state.E * 2 / (3 * params.N * params.k);
	}

	void print_xyz(std::ofstream& file) {
		file << params.N << std::endl;
		file << "Ar" << std::endl;
		for (int i = 0; i < params.N; i++)
			file << "Ar " << state.r[i].x() << " " << state.r[i].y() << " " << state.r[i].z() << std::endl;
	}

	void print_state(std::ofstream& file) {
		file << state.t << "," << state.E << "," << state.V << "," << state.H << "," << state.T << "," << state.P << std::endl;
	}

private:
	Parameters params;
	State state;
};



int main(int argc, const char* argv[]){
	time_t start,end;
	time(&start);
	Simulation sim(argv[1]);
	sim.run();
	time(&end);
	std::cout << "Time: " << difftime(end, start) << " s" << std::endl;
	return 0;
}
