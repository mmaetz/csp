#include <vector>
#include <random>
#include "functions.h"
#include "defn.h"
#include <iostream>
#include <cmath>

void initalize(std::vector< std::vector< std::vector<char> > >& cluster, const int seed, const double p, const double emax, double& ed, const int J)
{
	int icounter;
	int jcounter;
	int kcounter;

	std::uniform_int_distribution<> pos(0,N-1);
	std::mt19937 gen_pos(seed);

  for (icounter=0; icounter<N; icounter++)
	{
		for (jcounter=0; jcounter<N; jcounter++) 
		{
			for (kcounter=0; kcounter<N; kcounter++) 
			{
				cluster[icounter][jcounter][kcounter]=1;
			}
		}
	}
	int i,j,k;
	double this_energy;
	while(emax > ed + this_energy)
	{
		i = pos(gen_pos);
		j = pos(gen_pos);
		k = pos(gen_pos);
		this_energy = energy_diff(cluster, i, j, k, J);
		if( this_energy > 0 )
		{
			flip(cluster, i, j, k);
			ed += this_energy;
		}
	}
}

int energy_diff(std::vector< std::vector< std::vector<char> > >& cluster, int i, int j, int k, const int J)
{
	double energy = 0, energy1 = 0, energy2 = 0;
	energy1 += J*cluster[i][j][k]*cluster[period(i+1)][j][k];
	energy1 += J*cluster[i][j][k]*cluster[period(i-1)][j][k];
	energy1 += J*cluster[i][j][k]*cluster[i][period(j-1)][k];
	energy1 += J*cluster[i][j][k]*cluster[i][period(j+1)][k];
	energy1 += J*cluster[i][j][k]*cluster[i][j][period(k+1)];
	energy1 += J*cluster[i][j][k]*cluster[i][j][period(k-1)];

	energy2 -= J*cluster[i][j][k]*cluster[period(i+1)][j][k];
	energy2 -= J*cluster[i][j][k]*cluster[period(i-1)][j][k];
	energy2 -= J*cluster[i][j][k]*cluster[i][period(j-1)][k];
	energy2 -= J*cluster[i][j][k]*cluster[i][period(j+1)][k];
	energy2 -= J*cluster[i][j][k]*cluster[i][j][period(k+1)];
	energy2 -= J*cluster[i][j][k]*cluster[i][j][period(k-1)];

	energy = energy1-energy2;
	return energy;
}

void flip(std::vector< std::vector< std::vector<char> > >& cluster, int i, int j, int k)
{
	cluster[i][j][k] = -cluster[i][j][k];
}

int energy_tot(std::vector< std::vector< std::vector<char> > >& cluster, const int J)
{
	int i, j, k;
	double energy;

	std::vector< std::vector< std::vector<bool> > > visited(N, std::vector< std::vector<bool> >(N, std::vector<bool>(N)));

  for (i=0; i<N; i++)
	{
		for (j=0; j<N; j++) 
		{
			for (k=0; k<N; k++) 
			{
				visited[i][j][k] = 1;
//				if(!visited[i][j][k])
//					cluster[i][j][k]=1;
				if(!visited[period(i+1)][j][k]);
					energy -= J*cluster[i][j][k]*cluster[period(i+1)][j][k];
				if(!visited[period(i-1)][j][k]);
					energy -= J*cluster[i][j][k]*cluster[period(i-1)][j][k];
				if(!visited[i][period(j-1)][k]);
					energy -= J*cluster[i][j][k]*cluster[i][period(j-1)][k];
				if(!visited[i][period(j+1)][k]);
					energy -= J*cluster[i][j][k]*cluster[i][period(j+1)][k];
				if(!visited[i][j][period(k+1)]);
					energy -= J*cluster[i][j][k]*cluster[i][j][period(k+1)];
				if(!visited[i][j][period(k-1)]);
					energy -= J*cluster[i][j][k]*cluster[i][j][period(k-1)];
			}
		}
	}

	return energy;
}

double magn_tot(std::vector< std::vector< std::vector<char> > >& cluster)
{
	int i, j, k;
	double magnetization = 0;
  for (i=0; i<N; i++)
	{
		for (j=0; j<N; j++) 
		{
			for (k=0; k<N; k++) 
			{
					magnetization += cluster[i][j][k];
			}
		}
	}
	return magnetization/double(N);
}


int period(int index)
{
	switch (index) {
		case -1: 
			return N-1;
		case N:
			return 0;
		default:
			return index;
	}
}

const std::vector<double> prob_list(double beta)
{
	std::vector<double> list;
//	list.push_back(-exp(4.0/beta));
	list.push_back( std::min(double(1), exp(-4.0/beta)) );
	list.push_back( std::min(double(1), exp(-8.0/beta)) );
	list.push_back( std::min(double(1), exp(-12.0/beta)) );
	return list;
}

double metropolis_prob(int energy_diff, const std::vector<double>& list)
{
	switch (energy_diff) {
		case 4:
			return list[0];
		case 8:
			return list[1];
		case 12:
			return list[2];
	}
}
