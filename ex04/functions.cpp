#include <vector>
#include <random>
#include "functions.h"
#include "defn.h"
#include <iostream>
#include <cmath>

void initalize(std::vector< std::vector< std::vector<char> > >& cluster, const int seed, const double p)
{
	int icounter;
	int jcounter;
	int kcounter;

	std::mt19937 gen(seed);
	std::uniform_real_distribution<double> dis(0,1);

  for (icounter=0; icounter<N; icounter++)
	{
		for (jcounter=0; jcounter<N; jcounter++) 
		{
			for (kcounter=0; kcounter<N; kcounter++) 
			{
				if(dis(gen) < p)
				{
					cluster[icounter][jcounter][kcounter]=1;
				}
				else
					cluster[icounter][jcounter][kcounter]=-1;
			}
		}
	}
}

int energy_diff(std::vector< std::vector< std::vector<char> > >& cluster, int i, int j, int k, int J)
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

int energy_tot(std::vector< std::vector< std::vector<char> > >& cluster, int J)
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

void wolff(std::vector< std::vector< std::vector<char> > >& cluster, int i, int j, int k, const int seed)
{
	std::vector< std::vector< std::vector<char> > > tovisit(N, std::vector< std::vector<char> >(N, std::vector<char>(N)));
	std::mt19937 gen_pos(seed_pos);
	std::uniform_int_distribution<> pos(0,N-1);

	i = pos(gen_pos);
	j = pos(gen_pos);
	k = pos(gen_pos);

}

std::vector< std::vector<char> > add_neighbor( std::vector< std::vector< std::vector<char> > >& tovisit, int i, int j, int k)
{
	std::vector<char> tempvec(3);
	std::vector< std::vector<char> > list;
	if(tovisit[period(i+1)][j][k] == 0)
	{
		tovisit[period(i+1)][j][k]++;
		tempvec[0] = i+1;
		tempvec[1] = j;
		tempvec[2] = k;
		list.push_back(tempvec);
	}
	if(tovisit[period(i-1)][j][k] == 0)
	{
		tovisit[period(i-1)][j][k]++;
		tempvec[0] = period(i-1);
		tempvec[1] = j;
		tempvec[2] = k;
		list.push_back(tempvec);
	}
	if(tovisit[i][period(j+1)][k] == 0)
	{
		tovisit[i][period(j+1)][k]++;
		tempvec[0] = i;
		tempvec[1] = period(j+1);
		tempvec[2] = k;
		list.push_back(tempvec);
	}
	if(tovisit[i][period(j-1)][k] == 0)
	{
		tovisit[i][period(j-1)][k]++;
		tempvec[0] = i;
		tempvec[1] = period(j-1);
		tempvec[2] = k;
		list.push_back(tempvec);
	}
	if(tovisit[i][j][period(k+1)] == 0)
	{
		tovisit[i][j][period(k+1)]++;
		tempvec[0] = i;
		tempvec[1] = j;
		tempvec[2] = period(k+1);
		list.push_back(tempvec);
	}
	if(tovisit[i][j][period(k-1)] == 0)
	{
		tovisit[i][j][period(k-1)]++;
		tempvec[0] = i;
		tempvec[1] = j;
		tempvec[2] = period(k-1);
		list.push_back(tempvec);
	}
	return list;
}
