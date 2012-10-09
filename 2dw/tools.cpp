#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
const double pi = 3.1415;
double defsource(double x, double y)
{
    return 0;
}
void **matrix(int row, int col, int num_bytes)
  {
  int      i, num;
  char     **pointer, *ptr;

  pointer = new char* [row];
  if(!pointer) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for "<< row << "row addresses !" << endl;
    return NULL;
  }
  i = (row * col * num_bytes)/sizeof(char);
  pointer[0] = new char [i];
  if(!pointer[0]) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for address to " << i << " characters !" << endl;
    return NULL;
  }
  ptr = pointer[0];
  num = col * num_bytes;
  for(i = 0; i < row; i++, ptr += num )   {
    pointer[i] = ptr; 
  }

  return  (void **)pointer;

}
void free_matrix(void **matr)
{
	delete[] (char*) matr[0];
	delete[] matr;
}
void print_matrix(int n, int m,double** matr, const char * file, int skip)
{
	ofstream ofs(file);
	for (int i=0;i<n;i++)
	{
		if (i % skip==0) {
        for (int j=0;j<m;j++) {
			if (j %skip == 0)
				ofs << matr[i][j] << " ";
		}
		}
		ofs << endl;
	}
	ofs.close();
}
std::string ZeroPadNumber(int num)
{
    std::ostringstream ss;
    ss << std::setw(7) << std::setfill( '0' ) << num;
    return ss.str();
}
