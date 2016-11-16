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

double sum_recovered(vector<vd> recovered_species, double b, double rho, double k, double h, double i) {
	for (int i = 0; i < k; i++){
		total += recovered_species[h][i] * b * (1 - exp(-rho));
	}
	return total;
}

double calculate_susceptible_population(double g, double b, double term1, double sh, vector<vd> recovered_species, double t, double rho, double k, int h, int i){
	return g + b * exp(-term1) * sh + sum_recovered(recovered_species, b, rho, k, h, i);
}

double calculate_infected_population(double sh, double b, double term1, double term2, vector<vd> infected_species, double t, double gama, int h, int i){
	return sh * b * (1 - exp(- term1)) * term2 + infected_species[h][i] * b * exp(-gama);
}

double calculate_recovered(double b, double rho, double gama, vector<vd> infected_species, vector<vd> recovered_species) {
	return recovered_species * b * exp(-rho) + infected_species * b * (1 - exp(-gama));
}

int main() {
	double m = 0.03, alfa = 0.2, rho = 0.3;
	int h = 0, k = 0, t = 0;
	
	cout << "¿Cuál es el número de especies? " << endl;
	cin >> h;
	cout << "¿Cuál es el número de virus? " << endl;
	cin >> k;
	cout << "¿Cuantos días investigo? " << endl;
	cin >> t;

	vd s = vd(h, 1);
	vvd infected = vvd(h, vd(k, 1.0));
	vvd recovered = vvd(h, vd(k, 0.0));
	vvd beta = vvd(h, vd(k, 0.5));
	vd n = vd(h, 1);

	cout << "Escribe el total de cada población" << endl;

	for (int i = 0; i < h; i++) {
		cin >> n[i];
	}

	for (int j = 0; j < k; j++) {

	}

	return 0;
}