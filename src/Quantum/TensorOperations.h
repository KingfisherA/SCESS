#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

#define Slice2 slice<Eigen::array<Eigen::Index, 2>, Eigen::array<Eigen::Index, 2>>
#define Slice4 slice<Eigen::array<Eigen::Index, 4>, Eigen::array<Eigen::Index, 4>>

#define Shuffle2 shuffle<Eigen::array<Eigen::Index, 2>>
#define Shuffle4 shuffle<Eigen::array<Eigen::Index, 4>>

#define Contract1(other, dimensions) contract(other, Eigen::array<Eigen::IndexPair<Eigen::Index>,1>dimensions)
#define Contract2(other, dimensions) contract(other, Eigen::array<Eigen::IndexPair<Eigen::Index>,2>dimensions)
#define Contract3(other, dimensions) contract(other, Eigen::array<Eigen::IndexPair<Eigen::Index>,2>dimensions)
#define Contract4(other, dimensions) contract(other, Eigen::array<Eigen::IndexPair<Eigen::Index>,4>dimensions)

#define Reduce1 sum<Eigen::array<Eigen::Index, 1>>
#define Reduce2 sum<Eigen::array<Eigen::Index, 2>>
#define Reduce3 sum<Eigen::array<Eigen::Index, 3>>
#define Reduce4 sum<Eigen::array<Eigen::Index, 4>>