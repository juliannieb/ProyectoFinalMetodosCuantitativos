#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

typedef vector<double> vd;
typedef vector<vd> vvd;

double calculate_b(double m) {
	return 1 - m;
}

double calculate_g(double m, vector<double> k1, vector<double> n, int h) {
	return 2 * calculate_b(m) * n[h] * (1 - ((calculate_b(m) * n[h]) / k1[h]));
}

double calculate_term1(vector<vd> beta, vector<vd> infected_species, int h, int k){
	double sum = 0.0;
	for (int i = 0; i < k; i++){
		sum += (beta[h][i] * infected_species[h][i]);
	}
	return sum;
}

double calculate_term2(vector<vd> beta, vector<vd> infected_species, int h, int i, int k){
	return ((beta[h][i] * infected_species[h][i]) / calculate_term1(beta, infected_species, h, k));
}

double sum_recovered(vector<vd> recovered_species, double m, vvd rho, int h, double k) {
	int total = 0;
	for (int i = 0; i < k; i++){
		total += recovered_species[h][i] * calculate_b(m) * (1 - exp(-rho[h][i]));
	}
	return total;
}

double calculate_susceptible_population(vector<double> s, vector<vd> beta, vector<vd> recovered_species, vector<vd> infected_species, vvd rho, int k, int h, double m, vector<double> k1, vector<double> n){
	return calculate_g(m, k1, n, h) + calculate_b(m) * exp(-calculate_term1(beta, infected_species, h, k)) * s[h] + sum_recovered(recovered_species, m, rho, h, k);
}

double calculate_infected_population(vector<double> s, vector<vd> beta, vector<vd> infected_species, int h, int i, int k, vvd gama, double m){
	return s[h] * calculate_b(m) * (1 - exp(-calculate_term1(beta, infected_species, h, k))) * calculate_term2(beta, infected_species, h, i, k) + infected_species[h][i] * calculate_b(m) * exp(-gama[h][i]);
}

double calculate_recovered_population(vvd rho, vvd gama, vector<vd> infected_species, vector<vd> recovered_species, int h, int i, double m) {
	return recovered_species[h][i] * calculate_b(m) * exp(-rho[h][i]) + infected_species[h][i] * calculate_b(m) * (1 - exp(-gama[h][i]));
}

void print_matrix(vector<vd> matrix, int h, int k) {
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < k; j++) {
			if (j == 0) {
				printf("| %lf |", matrix[i][j]);
			} else {
				printf(" %lf |", matrix[i][j]);
			}
		}
		printf("\n");
	}
}

void print_vector(vector<double> vector_aux, int h) {
	for (int i = 0; i < h; i++) {
		if (i == 0) {
			printf("| %lf |", vector_aux[i]);
		} else {
			printf(" %lf |", vector_aux[i]);
		}
	}
}

void fill_random_matrix(vvd &matrix) {
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++) {
			double random_val = (double)(rand() % 500) / 1000;
			matrix[i][j] = random_val;
		}
	}
}

void fill_matrix(vvd &matrix){
	double val = 0;
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++) {
			printf("| mat[%i][%i]: |", i,j);
			cin >> val;
			matrix[i][j] = val;
		}
	}
}

int main() {
	srand (time(NULL));
	double m = 0.03, alfa = 0.2;
	int h = 0, k = 0, t = 0, species = 0;
	
	cout << "\nWhat's the amount of species? " << endl;
	cin >> h;
	cout << "\nWhat's the amount of virus? " << endl;
	cin >> k;
	cout << "\nHow many days do I investigate? " << endl;
	cin >> t;

	species = h;

	vvd infected = vvd(h, vd(k, 1.0));
	vvd recovered = vvd(h, vd(k, 0.0));
	vvd beta = vvd(h, vd(k, 0.5));
	vvd rho = vvd(h, vd(k, 0.5));
	vvd gama = vvd(h, vd(k, 0.5));

	//fill_random_matrix(beta);
	//fill_random_matrix(rho);
	//fill_random_matrix(gama);
	cout << "\nEnter the values of 'beta' \n" << endl;
	fill_matrix(beta);
	cout << "\nEnter the values of 'rho' " << endl;
	fill_matrix(rho);
	cout << "\nEnter the values of 'gama' " << endl;
	fill_matrix(gama);

	vd s = vd(h, 1);
	vd n = vd(h, 1);
	vd k1 = vd(h, 1);

	cout << "\nEnter the total of each population" << endl;

	for (int h = 0; h < species; h++) {
		cin >> n[h];
		k1[h] = n[h];
		s[h] = n[h];
	}

	printf("\nDay 0\n");

	printf("Susceptible table\n");
	print_vector(s, species);
	printf("\n\n");

	printf("Infected table\n");
	print_matrix(infected, species, k);
	printf("\n");

	printf("Recovered table\n");
	print_matrix(recovered, species, k);
	printf("\n");

	for (int a = 1; a <= t; a++) {

		printf("Day %d\n", a);

		for (int h = 0; h < species; h++) {
			s[h] = calculate_susceptible_population(s, beta, recovered, infected, rho, k, h, m, k1, n);

			for (int i = 0; i < k; i++) {
				infected[h][i] = calculate_infected_population(s, beta, infected, h, i, k, gama, m);

				recovered[h][i] = calculate_recovered_population(rho, gama, infected, recovered, h, i, m);
			}
		}

		printf("Susceptible table\n");
		print_vector(s, species);
		printf("\n\n");

		printf("Infected table\n");
		print_matrix(infected, species, k);
		printf("\n");

		printf("Recovered table\n");
		print_matrix(recovered, species, k);
		printf("\n");
	}

	return 0;
}