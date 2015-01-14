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
	double beta = 5;

	const int seed_occ = 1;
	const int seed_pos = 2;
	const int seed_en = 3;
	const double p = 0.5;
	const double efac = 3;
	const double emax = efac*N*N*N;
	std::cout << "Emax: " << emax << std::endl;
	double ed = 0;

	std::mt19937 gen_pos(seed_pos);
	std::mt19937 gen_en(seed_en);
	std::uniform_real_distribution<double> en(0,1);
	std::uniform_int_distribution<> pos(0,N-1);

	std::vector< std::vector< std::vector<char> > > cluster(N, std::vector< std::vector<char> >(N, std::vector<char>(N)));


////	 magnetism 1 or -1
	const int J = 1;
	initalize(cluster, seed_occ, p, emax, ed, J);
	Print_lattice (cluster[0], N, N, ImageWidth, ImageHeight, "random1.ppm");
	
	const int steps = pow(10, 7);
	int i,j,k;


	double etot = energy_tot(cluster, J);
	std::cout << "Total energy: " << etot << std::endl;
	double mtot = magn_tot(cluster);
	std::cout << "Total magnetization: " << mtot << std::endl;
	double this_energy;

	const std::vector<double> metro_prob = prob_list(beta);

	for(icounter = 0; icounter < steps; icounter++)
	{
		i = pos(gen_pos);
		j = pos(gen_pos);
		k = pos(gen_pos);
		this_energy = energy_diff(cluster, i, j, k, J);
		if( this_energy < 0 )
		{
			flip(cluster, i, j, k);
			ed += this_energy;
		}
		else
		{
//			if( std::min(double(1),exp(-double(this_energy)/beta)) >= en(gen_en))
			if( emax >= ed - this_energy && ed - this_energy >= 0)
			{
				flip(cluster, i, j, k);
				ed += this_energy;
			}
		}

	}
	Print_lattice (cluster[0], N, N, ImageWidth, ImageHeight, "random2.ppm");
	etot = energy_tot(cluster, J);
	std::cout << "Total energy: " << etot << std::endl;
	mtot = magn_tot(cluster);
	std::cout << "Total magnetization: " << mtot << std::endl;
}
