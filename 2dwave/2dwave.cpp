#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
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

  } // end: function void **matrix()

void print_matrix(int n,double** A, const char * file)
{
	ofstream ofs(file);
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
			ofs << A[i][j] << " ";
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
class wavesolver {
public:
	wavesolver(double(*_c)(double x, double y), double (*_initc)(double x,double y),double** _u, double _dt = 0.0001, double _h = 0.001, double _T=1, double _Lx = 1, double _Ly = 1)
	{
		initc = _initc;
		c = _c; 
		Lx = _Lx; Ly = _Ly; dt = _dt; h = _h; T=_T;
		Nx = (int) ceil(Lx/h);
		Nt = (int) ceil(T/dt);
		u = _u;
		current_time_step = 0;
	}

	int nextstep() {
		if (current_time_step==0) {
			init();
			return 1;
		}
		for (int i=1; i<Nx-1; i++) {
			for (int j=1; j<Nx-1; j++) {
				unext[i][j] = 2*u[i][j] - ulast[i][j] + (dt*dt)/(h*h)*(u[i+1][j] - 4*u[i][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]);
			}
		}
		for (int i=1; i<Nx-1; i++) {
			for (int j=1; j<Nx-1; j++) { 
				ulast[i][j] = u[i][j];
				u[i][j] = unext[i][j];
			}
		}
		return ++current_time_step;
	}
	double getTstep() {
		return Nt;
	}
private:
	double T, Lx, Ly, dt,h;
	int Nx,Ny,Nt;
	double **u;
	double **unext;
	double **ulast;
	int current_time_step;
	double(*c)(double x, double y);
	double(*initc)(double x, double y);
	void init() {
		if (current_time_step!=0)
			return;
		unext = (double**)matrix(Nx,Nx,sizeof(double));
		ulast = (double**)matrix(Nx,Nx,sizeof(double));
		for (int i=0; i<Nx; i++) {
			for (int j=0; j<Nx; j++) {
				u[i][j] = 0;
				unext[i][j] = 0;
				ulast[i][j] = (*initc)(h*i,h*j);
			}
		}
		// first timestep n = 0->1
		double temp;
		for (int i=1; i<Nx-1; i++) {
			for (int j=1; j<Nx-1; j++) {
				temp = ulast[i][j] + (dt)/(2*h*h)*(ulast[i+1][j] - 4*ulast[i][j] + ulast[i-1][j] + ulast[i][j+1] + ulast[i][j-1]);
				u[i][j] = 2*ulast[i][j] - temp + (dt*dt)/(h*h)*(ulast[i+1][j] - 4*ulast[i][j] + ulast[i-1][j] + ulast[i][j+1] + ulast[i][j-1]);
			}
		}
		current_time_step=1;

	}
};
const double pi = 3.1415;
double initcond(double x,double y)
{
	//return sin(pi*x)*sin(2*pi*y);
	return exp( -((x-0.5)*(x-0.5)/2 + (y-0.5)*(y-0.5))/2 );
}
double c(double x, double y)
{
	return 1;
}
int main()
{
	stringstream str;
	int write_delay = 10;
	double T = 1;
	double Lx = 1;
	double dt = 0.001;
	double h = 0.01;
	int Nx = (int) ceil(Lx/h);
	int Nt = (int) ceil(T/dt);
	double **u = (double**)matrix(Nx,Nx,sizeof(double));
	wavesolver wave1(&c,&initcond,u, dt, h, T, Lx);
	for (int n=0; n<wave1.getTstep(); n++) {
		wave1.nextstep();
		if (n % write_delay == 0) {
			str.str(std::string(""));
			str << "test.d" << ZeroPadNumber(n/write_delay);
			print_matrix(Nx,u,str.str().c_str());
		}
	}	
	/*	

	double **unext = (double**)matrix(Nx,Nx,sizeof(double));

	double **ulast = (double**)matrix(Nx,Nx,sizeof(double));

	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Nx; j++) {
			u[i][j] = 0;
			unext[i][j] = 0;
			ulast[i][j] = initcond(h*i,h*j);
		}
	}
	// first timestep n = 0
	double temp;
	for (int i=1; i<Nx-1; i++) {
		for (int j=1; j<Nx-1; j++) {
			temp = ulast[i][j] + (dt)/(2*h*h)*(ulast[i+1][j] - 4*ulast[i][j] + ulast[i-1][j] + ulast[i][j+1] + ulast[i][j-1]);
			u[i][j] = 2*ulast[i][j] - temp + (dt*dt)/(h*h)*(ulast[i+1][j] - 4*ulast[i][j] + ulast[i-1][j] + ulast[i][j+1] + ulast[i][j-1]);
		}
	}
	for (int n=1; n<Nt; n++) {
		for (int i=1; i<Nx-1; i++) {
			for (int j=1; j<Nx-1; j++) {
				unext[i][j] = 2*u[i][j] - ulast[i][j] + (dt*dt)/(h*h)*(u[i+1][j] - 4*u[i][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]);
				}
		}
		for (int i=1; i<Nx-1; i++) {
			for (int j=1; j<Nx-1; j++) { 
				ulast[i][j] = u[i][j];
				u[i][j] = unext[i][j];
			}
		}
		if (n % 100 == 0) {
		str.str(std::string(""));
		str << "test.d" << ZeroPadNumber(n/100);
		print_matrix(Nx,u,str.str().c_str());
		}
	} // time
	*/
	return 0;
}
