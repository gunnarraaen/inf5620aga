#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
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

void print_matrix(int n,double** matr, const char * file)
{
	ofstream ofs(file);
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
			ofs << matr[i][j] << " ";
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