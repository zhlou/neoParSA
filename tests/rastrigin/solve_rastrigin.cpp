#include <iostream>
#include "rastrigin.h"
#include "annealer.h"

int main(int argc, char **argv)
{
	rastrigin rst_problem(2);
	annealer rst_anneal(&rst_problem);

	return 0;
}
