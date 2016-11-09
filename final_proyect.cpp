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

int main() {
	cout << "Hola" << endl;
	return 0;
}