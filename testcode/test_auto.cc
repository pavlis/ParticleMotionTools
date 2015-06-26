#include <iostream>
#include "ParticleMotionError.h"
using namespace std;
/* test routine including some things I want to test*/
int main(int argc, char **argv)
{
    cout << "Test of C++ auto generation using ParticleMotionError"<<endl;
    cout << "Calling default constructor"<<endl;
    ParticleMotionError pm1;
    cout << "Default constuctor data:"<<endl
        << pm1 <<endl;
    cout << "Adding some data."<<endl<<"Result:"<<endl;
    pm1.dtheta_major=1.0;
    pm1.dmajornrm=3.3;
    pm1.ndgf_minor=4;
    cout << pm1<<endl;
    cout << "Calling copy constructor"<<endl;
    ParticleMotionError pm2(pm1);
    cout << "Result:"<<endl<<pm2<<endl;
    cout << "Calling automatically generated operator="<<endl;
    ParticleMotionError pm3;
    pm3=pm2;
    cout << "Result:"<<endl<<pm3<<endl;
}

