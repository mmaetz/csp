#include <iostream>
#include "latticeview.h"
#include "functions.h"
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>

main () {

	int icounter;
	int jcounter;

//	double kb = 8.6173324e-5;
//	double T = 300;
	double beta = 1e2;

	const int seed_occ = 1;
	const int seed_pos = 2;
	const int seed_en = 3;
	const double p = 0.5;

	std::mt19937 gen_pos(seed_pos);
	std::mt19937 gen_en(seed_en);
	std::uniform_real_distribution<double> en(0,1);
	std::uniform_int_distribution<> pos(0,N-1);

	std::vector< std::vector< std::vector<char> > > cluster(N, std::vector< std::vector<char> >(N, std::vector<char>(N)));


	initalize(cluster, seed_occ, p);
	Print_lattice (cluster[0], N, N, ImageWidth, ImageHeight, "random1.ppm");
	
	const int steps = pow(10, 6);
	int i,j,k;

////	 magnetism 1 or -1
	int J = 1;

	double etot = energy_tot(cluster, J);
	std::cout << etot << std::endl;
	double this_energy;

	for(icounter = 0; icounter < steps; icounter++)
	{
		i = pos(gen_pos);
		j = pos(gen_pos);
		k = pos(gen_pos);
		this_energy = energy_diff(cluster, i, j, k, J);
		if( this_energy <= 0 )
			flip(cluster, i, j, k);
		else
		{
			if( std::min(double(1),exp(-double(this_energy)/beta)) >= en(gen_en))
			{
				flip(cluster, i, j, k);
			}
		}

	}
	Print_lattice (cluster[5], N, N, ImageWidth, ImageHeight, "random2.ppm");
	etot = energy_tot(cluster, J);
	std::cout << etot << std::endl;
}
