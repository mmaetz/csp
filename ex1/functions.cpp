#include <vector>
#include <random>
#include "functions.h"
#include "defn.h"
#include <iostream>

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
				if(!visited[i][j][k])
					cluster[i][j][k]=1;
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

