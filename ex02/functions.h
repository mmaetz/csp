void initalize(std::vector< std::vector< std::vector<char> > >& cluster, const int seed, const double p);
int period(int index);
int energy_diff(std::vector< std::vector< std::vector<char> > >& cluster, int i, int j, int k, int J);
void flip(std::vector< std::vector< std::vector<char> > >& cluster, int i, int j, int k);
int energy_tot(std::vector< std::vector< std::vector<char> > >& cluster, int J);
const std::vector<double> prob_list(double beta);
double metropolis_prob(int energy_diff, const std::vector<double>& list);
double magn_tot(std::vector< std::vector< std::vector<char> > >& cluster);
