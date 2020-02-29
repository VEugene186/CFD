#ifndef PROFILE_H
#define PROFILE_H

#include "constants.h"

void setProfile(int nX, int nY, int & nXMin, int & nYMin, int & nXMax, int & nYMax, double ** x, double ** y, 
                double & xc, double & yc, int & iRib, int & jRib, const char * fileName);

#endif
