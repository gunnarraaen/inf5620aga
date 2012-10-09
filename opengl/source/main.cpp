#define ARMA_NO_DEBUG
#include <CWave.h>

int main(int a, char** b) {
    CWave p;

    // Wrapper: loads ini file and initialized all
    try {
        p.Initialize("wave.ini");
        p.Loop(p.Display, p.Update, std::string("MyBody"), a, b);
    }
    catch (string l) {
        cout << endl <<l << endl;
        return 1;
    }
    catch (...) {
        cout << "Unknown, fucked-up error.. should not happen!" << endl;
        return 1;
    }
    return 0;
}
