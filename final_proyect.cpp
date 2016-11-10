#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

typedef <vector<double> vd;
typedef vector<vd> vvd;

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

double calculate_susceptible_population(double s0, double m, double t, double k, double nh, vector<double> infected_species){
	if (t == 0) {
		return s0;
	}
	return calculate_g(m, t, k, nh) + calculate_b(m) * exp(-calculate_term1(calculate_b(m), infected_species, t)) * calculate_susceptible_population(s0, m, t-1, k, nh, infected_species) + ;
}

double calculate_infected_population(double calculate_susceptible_population, double s0, double m, double t, double k, double nh, vector<double> infected_species, ){

}

int main() {
	double m = 0.03, alfa = 0.2, rho = 0.3;
	int h = 0, k = 0, t = 0,;
	
	cout << "¿Cuál es el número de especies? " << endl;
	cin >> h;
	cout << "¿Cuál es el número de virus? " << endl;
	cin >> k;
	cout << "¿Cuantos días investigo? " << endl;
	cin >> t;

	vvd infected = vvd(h, vd(k, 1.0));
	vvd recovered = vvd(h, vd(k, 0.0));
	vvd beta = vvd(h, vd(k, 0.5));
	vd n = vd(h, 1);

	cout << "Escribe el total de cada población" << endl;
	for (int i = 0; i < h; i++) {
		cin >> n[i];
	}

	return 0;
}