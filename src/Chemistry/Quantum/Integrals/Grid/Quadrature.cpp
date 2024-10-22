//Derived from DFTCXX
//https://github.com/ifilot/dftcxx.git

#include<pch.h>
#include "Quadrature.h"

namespace SCESS 
{
	constexpr unsigned int Quadrature::nLebedevPoints[11];
	constexpr double Quadrature::lebedevCoefficients[][4];
}