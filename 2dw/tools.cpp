#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
using namespace boost::numeric::ublas;
const double pi = acos(-1);
// ----------------------------------
// print matrix to file
// ----------------------------------
void print_matrix(int n, int m, matrix<double> *matr, const char * file, int skip)
{
	ofstream ofs(file);
	for (int i=0;i<n;i++)
	{
		if (i % skip==0) {
        for (int j=0;j<m;j++) {
			if (j %skip == 0)
                ofs << (*matr)(i,j) << " ";
		}
		}
		ofs << endl;
	}
	ofs.close();
}
// ----------------------------------
// pad number with zeros
// ----------------------------------
std::string padnumber(int num,char pad)
{
    std::ostringstream ss;
    ss << std::setw(7) << std::setfill(pad) << num;
    return ss.str();
}
