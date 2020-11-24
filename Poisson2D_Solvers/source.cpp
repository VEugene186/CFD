#include "source.h"
#include <math.h>

double heatSource(double x, double y) {
    double dx = x - 0.9;
    double dy = y - 0.9;
    return 2000.0 * exp(-20.0 * (dx * dx + dy * dy));
}
