#include <iostream>
#include "rastrigin.h"
#include "annealer.h"
using namespace std;
int main(int argc, char **argv)
{
	rastrigin rst_problem(2);
	annealer rst_anneal(&rst_problem);
	cout << "The fininal energy is " << rst_anneal.loop() << endl;
	cout << "The solution is " << endl;
	rst_problem.print_solution(cout);

	return 0;
}
