#include <iostream>
#include "ParticleMotionError.h"
using namespace std;
void ParticleMotionError::zero()
{
    dtheta_major=0.0;  dphi_major=0.0;
    dtheta_minor=0.0;  dphi_minor=0.0;
    dmajornrm=0.0;  dminornrm=0.0;
    delta_rect=0.0;
    ndgf_major=0;
    ndgf_minor=0;
    ndgf_rect=0;
    ndgf_major_amp=0;
    ndgf_minor_amp=0;
}
ParticleMotionError::ParticleMotionError()
{
    this->zero();
}
ostream& operator<<(ostream& os, ParticleMotionError& pme)
{
    os << pme.dtheta_major<<" "
        <<pme.dphi_major<<" "
        <<pme.dtheta_minor<<" "
        <<pme.dphi_minor<<" "
        <<pme.dmajornrm<<" "
        <<pme.dminornrm<<" "
        <<pme.ndgf_major<<" "
        <<pme.ndgf_minor<<" "
        <<pme.ndgf_rect<<" "
        <<pme.ndgf_major_amp<<" "
        <<pme.ndgf_minor_amp;
    return(os);
}
