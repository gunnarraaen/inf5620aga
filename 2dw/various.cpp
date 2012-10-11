#include "tools.cpp"
#include "bitmap/EasyBMP.cpp"
// -------------------------------------------------------------
// Contains various initial conditions and wave-velocity coeffs
// -------------------------------------------------------------
// Sources
double no_source(double x, double y, double t)
{
    return 0;
}
double small_source(double x, double y, double t)
{
    return 2*exp(-((x-2)*(x-2) + (y-2)*(y-2))*180 );
}
// Initial condition
double initcond_exp(double x,double y)
{
    return 0;
    return exp( -((x-1)*(x-1) + (y-1)*(y-1))*70 );
}
double initcond_sin(double x, double y)
{
    return 2*sin(pi/4*x)*sin(pi/4*y);
}
double initcond_sin_exact(double x, double y, double t, double c)
{
    return 2*sin(pi/4.*x)*sin(pi/4.*y)*cos(pi/2.*c*c*t);
}
// Coeffs.
double coeff(double x, double y)
{
    if (x>1.6 && x<1.7)
    {
        if ((y>1.5 && y<1.6) || (y>2.4 && y<2.5))
            return 1;
        else
            return 0;
    }
    return 1;
}
double coeff2(double x, double y)
{
    return 1 - exp( - ((x-2)*(x-2) + (y-2)*(y-2))*20 ) - exp( -((x-3)*(x-3) + (y-3)*(y-3))*20 );
}
double coeff3(double x, double y)
{
    if (x<0 || y<0) return 0;
    int Lx = 4;
    int Ly = 4;
    static bool used_once = false;
    static BMP input;
    if (!used_once)
    {   // open image
        input.ReadFromFile("test4.bmp");
        used_once = true;
    }
    int new_x = x*input.TellHeight()/Lx;
    int new_y = y*input.TellWidth()/Ly;

    return 0.213*(int)input(new_y,new_x)->Red + 0.715*(int)input(new_y,new_x)->Green + 0.072*(int)input(new_y,new_x)->Blue;
}
