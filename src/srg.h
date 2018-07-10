#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;



static const int Bmax = 32;
static const vec B = {1.0, -1/2.0, 1/6.0, 0.0, -1/30.0, 
	                  0.0, 1/42.0, 0.0, -1/30.0, 0.0, 
	                  5/66.0, 0.0, -691/2730.0, 0.0, 7/6.0, 
	                  0.0, -3617/510.0, 0.0, 43867/798.0, 0.0, 
	                  -174611/330, 0.0, 854513/138.0, 0.0, -236364091/2730.0, 
	                  0.0, 8553103/6.0, 0.0, -23749461029/870, 0.0,
	                  8615841276005/14322, 0.0};