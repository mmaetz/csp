#include<iostream>
#include<fstream>
#include <vector>

#define N 100             //Lateral number of cells
#define ImageWidth 1000  //image width
#define ImageHeight 1000 //image height

using namespace std;

void Print_lattice (std::vector< std::vector<bool> >& vlat, const int &vlx, const int &vly, const int &vwidth, const int &vheight, const char* vfilename="output.ppm")
{
  int  i, j, k, l;
  int vw= vwidth/vlx, vh=vheight/vly;
  int r[2], g[2], b[2];

  r[0]= 255; g[0]= 255; b[0]= 255; //white
  r[1]=   0; g[1]=   0; b[1]=   0; //black

  ofstream out (vfilename);

  out << "P3" << endl;
  out << vw*vlx << " " << vh*vly << endl;
  out << "255" << endl;

  for (i=vly-1; i>=0; i--)
    for (j=0; j<vh; j++)
      for (k=0; k<vlx; k++)
      {
        for (l=0; l<vw; l++)
        { out << r[vlat[i][k]] << " " << g[vlat[i][k]] << " " << b[vlat[i][k]] << " ";
        }
      } 
      out << endl;

  out.close ();
}
