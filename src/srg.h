#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;



static const int Bmax = 22;
static const vec B = {1.0, -1/2.0, 1/6.0, 0.0, -1/30.0, 
	                  0.0, 1/42.0, 0.0, -1/30.0, 0.0, 5/66.0, 
	                  0.0, -691/2730.0, 0.0, 7/6.0, 0.0, -3617/510.0, 
	                  0.0, 43867/798.0, 0.0, -174611/330, 0.0}


