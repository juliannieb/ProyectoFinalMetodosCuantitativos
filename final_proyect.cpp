#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

double calculate_b(double m) {
	return 1 - m;
}

double calculate_g(double m, double t, double k, double nh) {
	double g = 2 * calculate_b(m) * nh * t * (1 - ((calculate_b(m) * nh * t) / k));
	return g;
}

double calculate_term1(double b, vector<double> infected_species, double t){
	double sum = 0.0;
	int k = infected_species.size();
	for (int i = 0; i < k; i++){
		sum += (b * infected_species[i] * t);
	}
	return sum;
}

double calculate_term2(double b, vector<double> infected_species, double t){
	double sum = ((b * infected_species[0] * t)/calculate_term1(b, infected_species, t));
	return sum;
}

int main() {
	cout << "Hola" << endl;
	return 0;
}