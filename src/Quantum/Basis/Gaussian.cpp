#include<pch.h>
#include "Gaussian.h"

//
//
const std::unordered_map<SCESS::OrbitalType, Eigen::Vector3i> SCESS::GaussianCollection::sAngMomMap = {
            {OrbitalType::S,    {0, 0, 0}} ,
            //
            {OrbitalType::P_X , {1, 0, 0} } ,
            {OrbitalType::P_Y , {0, 1, 0} } ,
            {OrbitalType::P_Z , {0, 0, 1} } ,
            //
            {OrbitalType::D_X2 , {2, 0, 0} },
            {OrbitalType::D_Y2 , {0, 2, 0} },
            {OrbitalType::D_Z2 , {0, 0, 2} },
            {OrbitalType::D_XY , {1, 1, 0} },
            {OrbitalType::D_XZ , {1, 0, 1} },
            {OrbitalType::D_YZ , {0, 1, 1} },
            //
            {OrbitalType::F_X3 , {3, 0, 0} },
            {OrbitalType::F_Y3 , {0, 3, 0} },
            {OrbitalType::F_Z3 , {0, 0, 3} },
            {OrbitalType::F_X2_Y , {2, 1, 0} },
            {OrbitalType::F_X_Y2 , {1, 2, 0} },
            {OrbitalType::F_X2_Z , {2, 0, 1} },
            {OrbitalType::F_X_Z2 , {1, 0, 2} },
            {OrbitalType::F_Y2_Z , {0, 2, 1} },
            {OrbitalType::F_Y_Z2 , {0, 1, 2} },
            {OrbitalType::F_XYZ , {1, 1, 1} }
};
