#include "tools.cpp"
#include "bitmap/EasyBMP.cpp"
// -------------------------------------------------------------
// Contains various initial conditions, velocity fields and sources
// -------------------------------------------------------------
// Sources
double no_source(double x, double y, double t) {
    return 0;
}
double small_source(double x, double y, double t) {
    return exp(-t)*exp(-((x-2)*(x-2) + (y-2)*(y-2))*180 );
}
// Initial conditions
double initcond_exp(double x,double y) {
    return exp( -((x-1)*(x-1) + (y-1)*(y-1))*70 );
}
// Vel. fields
double coeff(double x, double y) {
    if (x>1.6 && x<1.7)
    {
        if ((y>1.5 && y<1.6) || (y>2.4 && y<2.5))
            return 1;
        else
            return 0;
    }
    return 1;
}
double coeff2(double x, double y) {
    return 1 - exp( - ((x-2)*(x-2) + (y-2)*(y-2))*20 );
}
double Vzero(double x, double y) {
    return 0;
}
double box_initc(double x, double y) {
    if ( 1.5 < x && x < 2.5)
        return 1;
    else
        return 0;
}
// -----------------------
// standing wave
// -----------------------
double omega = pi;
double b = 0.01;
double Lx = 4;
double Ly = 4;
int mx = 4;
int my = 4;
// Source term
double exact_source(double x, double y, double t) {
    return exp(-b*t)*cos(mx*pi*x/Lx)*cos(my*pi*y/Ly)*(b*omega*sin(omega*t) + (-omega*omega + (mx*pi/Lx)*(mx*pi/Lx) + (my*pi/Ly)*(my*pi/Ly))*cos(omega*t));
}
// V = du/dt at t=0
double exact_V(double x, double y) {
    return -b*cos(mx*pi*x/Lx)*cos(my*pi*y/Ly);
}
// Initial condition
double exact_initc(double x, double y) {
    return cos(mx*pi*x/Lx)*cos(my*pi*y/Ly);
}
// Velocity field
double exact_coeff(double x, double y) {
    return 1;
}
// Exact solution
double exact_sol(double x, double y, double t) {
    return exp(-b*t)*cos(mx*pi*x/Lx)*cos(my*pi*y/Ly)*cos(omega*t);
}
// -----------------------

// -----------------------
// vel. field from bitmap
// -----------------------
double coeff3(double x, double y) {
    if (x<0 || y<0) return 0;
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
