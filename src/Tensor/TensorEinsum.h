#pragma once

#include "TensorOperations.h"
#include <array>
#include <unsupported/Eigen/CXX11/Tensor>

namespace SCESS {
    //
    //
#pragma region Meta Classes

    template<class IndexI, class IndexJ>
    struct is_reduction_op {
        static constexpr bool value = IndexI::N != IndexJ::N;
    };

    template<class IndexI, class IndexJ>
    struct is_shuffle_op {
        static constexpr bool value = IndexI::N == IndexJ::N;
    };

#pragma endregion
    //
    //
#pragma region Shuffle

    template<Idx N>
    struct ShuffleIdxs {
        std::array<Idx, N> shuffle{};
        Idx nShuffle{0};
    };

    template<class Idxs>
    constexpr bool IsShuffleRedundant(Idxs idxs) {
        bool isShuffleRedundant = true;

        for (Idx i = 0; i < idxs.shuffle.size(); i++)
            if (i != idxs.shuffle[i]) {
                isShuffleRedundant = false;
                break;
            }

        return isShuffleRedundant;
    }

    template<class IndexI, class IndexJ>
    constexpr ShuffleIdxs<IndexI::N> GenerateShuffleIdxs() {
        static_assert(IndexI::N == IndexJ::N, "IndexI and IndexJ must have the same size");
        ShuffleIdxs<IndexI::N> idxs;

        for (Idx j = 0; j < IndexJ::values.size(); ++j)
            for (Idx i = 0; i < IndexI::values.size(); ++i)
                if (IndexI::values[i] == IndexJ::values[j]) {
                    idxs.shuffle[j] = i;
                    ++idxs.nShuffle;
                }


        return idxs;
    }

    //Einsum for a shuffle operator
    template<class IndexI, class IndexJ,
        class TensorType,
        std::enable_if_t<is_shuffle_op<IndexI, IndexJ>::value, bool>  = true>
    auto Einsum(TensorType t1) {
        static_assert(IndexI::N == IndexJ::N, "IndexI and IndexJ must have the same size");

        constexpr auto idxs = GenerateShuffleIdxs<IndexI, IndexJ>();
        static_assert(IndexI::N == idxs.nShuffle, "IndexI and IndexJ do not have the same shuffling indices");

        if constexpr (IsShuffleRedundant<std::remove_reference_t<decltype(idxs)> >(idxs))
            return std::ref(t1);
        else
            return t1.template shuffle<std::remove_reference_t<decltype(idxs.shuffle)> >(idxs.shuffle);
    }

#pragma endregion
    //
    //
#pragma region Reduction

    template<Idx N_I, Idx N_J>
    struct ReductionIdxsTmp {
        std::array<Idx, N_I> reduction{};
        std::array<Idx, N_J> shuffle{};
        Idx nReduce{0};
        Idx nShuffle{0};
    };

    template<Idx nReduce, Idx nShuffle>
    struct ReductionIdxs {
        std::array<Idx, nReduce> reduction{0};
        std::array<Idx, nShuffle> shuffle{0};
    };

    template<class IndexI, class IndexJ,
        class ReturnT = ReductionIdxsTmp<2 * IndexI::N, IndexJ::N> >
    constexpr ReturnT GenerateReductionIdxsTmp() {
        static_assert(IndexI::N > IndexJ::N, "IndexJ must be smaller than IndexI in order to perform a reduction");
        ReturnT idxs;
        std::array<Idx, 2 * IndexI::N> reductionIdxsTmp{};
        Idx p = 0;

        for (Idx i = 0; i < IndexI::N; ++i)
            for (Idx j = 0; j < IndexI::N; ++j)
                if (IndexI::values[i] == IndexI::values[j] && i != j) {
                    reductionIdxsTmp[p] = IndexI::values[i];
                    idxs.reduction[p++] = i;
                }

        idxs.nReduce = p;
        std::array<Idx, IndexI::N> newIndexI;
        p = 0;

        for (Idx i = 0; i < IndexI::N; ++i) {
            bool shouldWrite = true;
            for (Idx j = 0; j < idxs.nReduce; ++j)
                shouldWrite &= reductionIdxsTmp[j] != IndexI::values[i];

            if (shouldWrite)
                newIndexI[p++] = IndexI::values[i];
        }

        Idx j = 0;

        for (; j < IndexJ::N; ++j)
            for (Idx i = 0; i < p; ++i)
                if (IndexJ::values[j] == newIndexI[i]) {
                    idxs.shuffle[j] = i;
                    break;
                }

        idxs.nShuffle = j;

        return idxs;
    }

    template<class IndexI, class IndexJ,
        auto reduce = GenerateReductionIdxsTmp<IndexI, IndexJ>()>
    constexpr ReductionIdxs<reduce.nReduce, reduce.nShuffle> GenerateReductionIdxs() {
        ReductionIdxs<reduce.nReduce, reduce.nShuffle> idxs;

        for (Idx i = 0; i < reduce.nReduce; ++i)
            idxs.reduction[i] = reduce.reduction[i];

        for (Idx i = 0; i < reduce.nShuffle; ++i)
            idxs.shuffle[i] = reduce.shuffle[i];

        return idxs;
    }

    //Einsum for a reduction and shuffle operator
    template<class IndexI, class IndexJ,
        class TensorType,
        std::enable_if_t<is_reduction_op<IndexI, IndexJ>::value, bool>  = true>
    auto Einsum(TensorType t1) {
        constexpr auto idxs = GenerateReductionIdxs<IndexI, IndexJ>();

        if constexpr (IsShuffleRedundant<decltype(idxs)>(idxs))
            return t1.template sum<std::remove_reference_t<decltype(idxs.reduction)> >(idxs.reduction);
        else
            return t1.template sum<std::remove_reference_t<decltype(idxs.reduction)> >(idxs.reduction)
                    .template shuffle<std::remove_reference_t<decltype(idxs.shuffle)> >(idxs.shuffle);
    }

#pragma endregion
    //
    //
#pragma region Contraction

    template<Idx N_I, Idx N_K>
    struct ContractIdxsTmp {
        std::array<IdxP, N_I> contract{};
        std::array<Idx, N_K> shuffle{};
        Idx nContract{0};
        Idx nShuffle{0};
    };

    template<Idx nContract, Idx nShuffle>
    struct ContractIdxs {
        std::array<IdxP, nContract> contract{};
        std::array<Idx, nShuffle> shuffle{};
    };

    template<class IndexI, class IndexJ, class IndexK>
    constexpr ContractIdxsTmp<IndexI::N, IndexK::N> GenerateContractIdxsTmp() {
        ContractIdxsTmp<IndexI::N, IndexK::N> idxs;
        std::array<IdxP, IndexI::N> idxContractTmp{{{0, 0}}};
        Idx p = 0;

        for (Idx i = 0; i < IndexI::values.size(); ++i)
            for (Idx j = 0; j < IndexJ::values.size(); ++j)
                if (IndexI::values[i] == IndexJ::values[j]) {
                    idxContractTmp[p] = {IndexI::values[i], IndexJ::values[j]};
                    idxs.contract[p++] = {i, j};
                }

        idxs.nContract = p;

        std::array<Idx, IndexI::N + IndexJ::N> newIndexI;
        p = 0;

        for (Idx i = 0; i < IndexI::N; ++i) {
            bool shouldWrite = true;
            for (Idx j = 0; j < idxs.nContract; ++j)
                shouldWrite &= idxContractTmp[j].first != IndexI::values[i];

            if (shouldWrite)
                newIndexI[p++] = IndexI::values[i];
        }

        for (Idx i = 0; i < IndexJ::N; ++i) {
            bool shouldWrite = true;
            for (Idx j = 0; j < idxs.nContract; ++j)
                shouldWrite &= idxContractTmp[j].second != IndexJ::values[i];

            if (shouldWrite)
                newIndexI[p++] = IndexJ::values[i];
        }

        Idx k = 0;

        for (; k < IndexK::N; ++k)
            for (Idx i = 0; i < p; ++i)
                if (IndexK::values[k] == newIndexI[i]) {
                    idxs.shuffle[k] = i;
                    break;
                }

        idxs.nShuffle = k;

        return idxs;
    }

    template<class IndexI, class IndexJ, class IndexK,
        ContractIdxsTmp<IndexI::N, IndexK::N> contract = GenerateContractIdxsTmp<IndexI, IndexJ, IndexK>()>
    constexpr ContractIdxs<contract.nContract, contract.nShuffle> GenerateContractIdxs() {
        ContractIdxs<contract.nContract, contract.nShuffle> idxs;

        static_assert((IndexI::N + IndexJ::N - 2 * contract.nContract) == IndexK::N,
                      "Contraction dimensions and IndexK must be the same size");
        static_assert(contract.nShuffle == IndexK::N, "Shuffle dimensions and IndexK must be the same size");


        for (Idx i = 0; i < contract.nContract; ++i)
            idxs.contract[i] = contract.contract[i];

        for (Idx i = 0; i < contract.nShuffle; ++i)
            idxs.shuffle[i] = contract.shuffle[i];


        return idxs;
    }

    //Einsum for contraction and shuffle
    template<class IndexI, class IndexJ, class IndexK,
        class TensorType1, class TensorType2>
    auto Einsum(TensorType1 t1, TensorType2 t2) {
        constexpr auto idxs = GenerateContractIdxs<IndexI, IndexJ, IndexK>();

        if constexpr (IsShuffleRedundant<decltype(idxs)>(idxs))
            return t1.template contract<std::remove_reference_t<TensorType2>,
                std::remove_reference_t<decltype(idxs.contract)> >(t2, idxs.contract);
        else
            return t1.template contract<std::remove_reference_t<TensorType2>,
                        std::remove_reference_t<decltype(idxs.contract)> >(t2, idxs.contract)
                    .template shuffle<std::remove_reference_t<decltype(idxs.shuffle)> >(idxs.shuffle);
    }

#pragma endregion
    //
    //
#pragma region Specialisation
    template<class IndexI, class IndexJ,
        class Scalar, int tensorRank>
    auto Einsum(Eigen::Tensor<Scalar, tensorRank> &t1) {
        return Einsum<IndexI, IndexJ, decltype(t1)>(t1);
    }

    template<class IndexI, class IndexJ,
        class Scalar, int tensorRank>
    auto Einsum(const Eigen::Tensor<Scalar, tensorRank> &t1) {
        return Einsum<IndexI, IndexJ, decltype(t1)>(t1);
    }

    template<class IndexI, class IndexJ, class IndexK,
        class Scalar1, int tensorRank1, class TensorType2>
    auto Einsum(Eigen::Tensor<Scalar1, tensorRank1> &t1, TensorType2 t2) {
        return Einsum<IndexI, IndexJ, IndexK, decltype(t1), decltype(t2)>(t1, t2);
    }

    template<class IndexI, class IndexJ, class IndexK,
        class Scalar1, int tensorRank1, class TensorType2>
    auto Einsum(const Eigen::Tensor<Scalar1, tensorRank1> &t1, TensorType2 t2) {
        return Einsum<IndexI, IndexJ, IndexK, decltype(t1), decltype(t2)>(t1, t2);
    }

    template<class IndexI, class IndexJ, class IndexK,
        class TensorType1, class Scalar2, int tensorRank2>
    auto Einsum(TensorType1 t1, Eigen::Tensor<Scalar2, tensorRank2> &t2) {
        return Einsum<IndexI, IndexJ, IndexK, decltype(t1), decltype(t2)>(t1, t2);
    }

    template<class IndexI, class IndexJ, class IndexK,
        class TensorType1, class Scalar2, int tensorRank2>
    auto Einsum(TensorType1 t1, const Eigen::Tensor<Scalar2, tensorRank2> &t2) {
        return Einsum<IndexI, IndexJ, IndexK, decltype(t1), decltype(t2)>(t1, t2);
    }

    template<class IndexI, class IndexJ, class IndexK,
        class Scalar1, int tensorRank1, class Scalar2, int tensorRank2>
    auto Einsum(Eigen::Tensor<Scalar1, tensorRank1> &t1, Eigen::Tensor<Scalar2, tensorRank2> &t2) {
        return Einsum<IndexI, IndexJ, IndexK, decltype(t1), decltype(t2)>(t1, t2);
    }

    template<class IndexI, class IndexJ, class IndexK,
        class Scalar1, int tensorRank1, class Scalar2, int tensorRank2>
    auto Einsum(const Eigen::Tensor<Scalar1, tensorRank1> &t1, const Eigen::Tensor<Scalar2, tensorRank2> &t2) {
        return Einsum<IndexI, IndexJ, IndexK, decltype(t1), decltype(t2)>(t1, t2);
    }

#pragma endregion
    //
    //
}
