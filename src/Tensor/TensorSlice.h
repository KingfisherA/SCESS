#pragma once

#include "TensorOperations.h"
#include <unordered_map>

namespace SCESS {
    // index -> {startIndex, size}
    using TensorSliceIndices = std::unordered_map<int, std::pair<Idx, Idx> >;

    template<class IndexI, class TensorType>
    auto Slice(const TensorSliceIndices &idxs, TensorType t1) {
        std::array<Idx, IndexI::N> startIdxs;
        std::array<Idx, IndexI::N> sizeArr;

        for (Idx i = 0; i < IndexI::N; i++) {
            const auto &pair = idxs.at(IndexI::values[i]);
            startIdxs[i] = pair.first;
            sizeArr[i] = pair.second;
        }

        return t1.template slice<decltype(startIdxs), decltype(sizeArr)>(startIdxs, sizeArr);
    }

    //Specialisations
    template<class IndexI, class Scalar, int tensorRank>
    auto Slice(const TensorSliceIndices &idxs, Eigen::Tensor<Scalar, tensorRank> &t1) {
        return Slice<IndexI, decltype(t1)>(idxs, t1);
    }

    template<class IndexI, class Scalar, int tensorRank>
    auto Slice(const TensorSliceIndices &idxs, const Eigen::Tensor<Scalar, tensorRank> &t1) {
        return Slice<IndexI, decltype(t1)>(idxs, t1);
    }
}
