#pragma once

#include "TensorOperations.h"
#include <array>
#include <type_traits>
#include <unsupported/Eigen/CXX11/Tensor>

//Compile-time Einstein summation over Eigen tensor expressions. Each operand is labelled by an Indices<...> list:
//  Einsum<Indices<I, J>, Indices<J, I> >(t)                         ij->ji       (shuffle)
//  Einsum<Indices<I, J, K>, Indices<K, I> >(t)                      ijk->ki      (reduction then shuffle)
//  Einsum<Indices<I, J>, Indices<J, K>, Indices<I, K> >(t1, t2)     ij,jk->ik    (contraction then shuffle)
//Labels shared by both operands of a contraction are summed over. A label repeated within one operand (a trace or
//diagonal) and a label kept in the output from both operands (a Hadamard product) are not supported.
//
//The result is an unevaluated Eigen expression. Like Eigen itself, plain tensors are referenced rather than copied,
//so a tensor operand must outlive the result; a temporary tensor is only safe within the same full expression.
namespace SCESS {
    //
    //
#pragma region Index Arithmetic

    namespace EinsumDetail {
        template<std::size_t N>
        constexpr Idx Find(const std::array<Idx, N> &labels, Idx label) {
            for (Idx i = 0; i < static_cast<Idx>(N); ++i)
                if (labels[i] == label)
                    return i;
            return -1;
        }

        template<std::size_t N>
        constexpr bool HasRepeatedIndices(const std::array<Idx, N> &labels) {
            for (Idx i = 0; i < static_cast<Idx>(N); ++i)
                if (Find(labels, labels[i]) != i)
                    return true;
            return false;
        }

        //Number of labels of IndexA which also appear in IndexB
        template<class IndexA, class IndexB>
        constexpr std::size_t CountShared() {
            std::size_t n = 0;
            for (Idx label: IndexA::values)
                n += Find(IndexB::values, label) >= 0;
            return n;
        }

        //Labels (or their positions) of IndexA which do or do not also appear in IndexB, in IndexA order
        template<class IndexA, class IndexB, bool shared, bool positions = false>
        constexpr auto SelectIndices() {
            constexpr std::size_t nShared = CountShared<IndexA, IndexB>();
            std::array<Idx, shared ? nShared : IndexA::N - nShared> selected{};
            std::size_t p = 0;
            for (Idx a = 0; a < IndexA::N; ++a)
                if ((Find(IndexB::values, IndexA::values[a]) >= 0) == shared)
                    selected[p++] = positions ? a : IndexA::values[a];
            return selected;
        }

        template<std::size_t N, std::size_t M>
        constexpr std::array<Idx, N + M> Concat(const std::array<Idx, N> &first, const std::array<Idx, M> &second) {
            std::array<Idx, N + M> result{};
            for (std::size_t i = 0; i < N; ++i)
                result[i] = first[i];
            for (std::size_t i = 0; i < M; ++i)
                result[N + i] = second[i];
            return result;
        }

        //Eigen shuffle taking dimensions labelled `from` to the order `to`: output dimension k is input dimension
        //perm[k]. A label of `to` missing from `from` gives -1.
        template<std::size_t N, std::size_t M>
        constexpr std::array<Idx, M> Permutation(const std::array<Idx, N> &from, const std::array<Idx, M> &to) {
            std::array<Idx, M> perm{};
            for (std::size_t k = 0; k < M; ++k)
                perm[k] = Find(from, to[k]);
            return perm;
        }

        //True if `to` (free of repeated labels) is a rearrangement of exactly the labels of `from`
        template<std::size_t N, std::size_t M>
        constexpr bool IsPermutationOf(const std::array<Idx, N> &from, const std::array<Idx, M> &to) {
            if (N != M)
                return false;
            for (Idx p: Permutation(from, to))
                if (p < 0)
                    return false;
            return true;
        }

        template<std::size_t N>
        constexpr bool IsIdentity(const std::array<Idx, N> &perm) {
            for (std::size_t i = 0; i < N; ++i)
                if (perm[i] != static_cast<Idx>(i))
                    return false;
            return true;
        }

        //A shuffle which leaves the innermost (contiguous) dimension in place can copy whole contiguous runs
        template<std::size_t N>
        constexpr bool KeepsInnermostDim(const std::array<Idx, N> &perm, int layout) {
            if constexpr (N == 0)
                return true;
            else {
                const Idx inner = layout == Eigen::ColMajor ? 0 : static_cast<Idx>(N) - 1;
                return perm[inner] == inner;
            }
        }

        //Pairs {dimension of IndexA, dimension of IndexB} over the labels shared by both operands
        template<class IndexA, class IndexB>
        constexpr auto ContractionPairs() {
            std::array<IdxP, CountShared<IndexA, IndexB>()> pairs{};
            std::size_t p = 0;
            for (Idx a = 0; a < IndexA::N; ++a)
                if (const Idx b = Find(IndexB::values, IndexA::values[a]); b >= 0)
                    pairs[p++] = {a, b};
            return pairs;
        }

        template<class TensorType>
        constexpr int Rank = Eigen::internal::traits<std::remove_cvref_t<TensorType> >::NumDimensions;

        template<class TensorType>
        constexpr int Layout = Eigen::internal::traits<std::remove_cvref_t<TensorType> >::Layout;
    }

#pragma endregion
    //
    //
#pragma region Shuffle & Reduction

    //Labels of IndexI missing from IndexJ are summed over, the remaining dimensions are shuffled into IndexJ order.
    //When IndexJ is a rearrangement of IndexI this is a pure shuffle, and when it is the same order the operand is
    //returned untouched (by reference for an lvalue).
    template<class IndexI, class IndexJ, class TensorType>
    decltype(auto) Einsum(TensorType &&t1) {
        using namespace EinsumDetail;
        static_assert(Rank<TensorType> == IndexI::N, "The rank of the tensor does not match IndexI");
        static_assert(!HasRepeatedIndices(IndexI::values) && !HasRepeatedIndices(IndexJ::values),
                      "Repeated indices within one operand (traces) are not supported by Einsum");

        constexpr auto reduced = SelectIndices<IndexI, IndexJ, false, true>();
        constexpr auto kept = SelectIndices<IndexI, IndexJ, true>();
        static_assert(IsPermutationOf(kept, IndexJ::values), "Every index of IndexJ must appear in IndexI");
        constexpr auto perm = Permutation(kept, IndexJ::values);

        if constexpr (reduced.size() > 0) {
            if constexpr (IsIdentity(perm))
                return t1.sum(reduced);
            else
                return t1.sum(reduced).shuffle(perm);
        } else if constexpr (!IsIdentity(perm))
            return t1.shuffle(perm);
        else if constexpr (std::is_lvalue_reference_v<TensorType>)
            return t1;
        else
            return std::remove_cvref_t<TensorType>(std::move(t1));
    }

#pragma endregion
    //
    //
#pragma region Contraction

    //Labels shared by IndexI and IndexJ are contracted; IndexK orders the remaining labels. Eigen lays out a
    //contraction as the free dimensions of its left operand followed by those of its right operand, so the operands
    //are swapped whenever that removes the output shuffle, or failing that lets it keep the innermost dimension.
    template<class IndexI, class IndexJ, class IndexK, class TensorType1, class TensorType2>
    auto Einsum(TensorType1 &&t1, TensorType2 &&t2) {
        using namespace EinsumDetail;
        static_assert(Rank<TensorType1> == IndexI::N, "The rank of the first tensor does not match IndexI");
        static_assert(Rank<TensorType2> == IndexJ::N, "The rank of the second tensor does not match IndexJ");
        static_assert(Layout<TensorType1> == Layout<TensorType2>, "Both tensors must have the same storage order");
        static_assert(!HasRepeatedIndices(IndexI::values) && !HasRepeatedIndices(IndexJ::values) &&
                      !HasRepeatedIndices(IndexK::values),
                      "Repeated indices within one operand (traces) are not supported by Einsum");

        constexpr auto freeI = SelectIndices<IndexI, IndexJ, false>();
        constexpr auto freeJ = SelectIndices<IndexJ, IndexI, false>();
        static_assert(IsPermutationOf(Concat(freeI, freeJ), IndexK::values),
                      "IndexK must contain exactly the indices which are not contracted");

        constexpr auto perm = Permutation(Concat(freeI, freeJ), IndexK::values);
        constexpr auto permSwapped = Permutation(Concat(freeJ, freeI), IndexK::values);
        constexpr int layout = Layout<TensorType1>;
        constexpr bool swap = !IsIdentity(perm) &&
                              (IsIdentity(permSwapped) ||
                               (KeepsInnermostDim(permSwapped, layout) && !KeepsInnermostDim(perm, layout)));

        if constexpr (swap) {
            auto contraction = t2.contract(t1, ContractionPairs<IndexJ, IndexI>());
            if constexpr (IsIdentity(permSwapped))
                return contraction;
            else
                return contraction.shuffle(permSwapped);
        } else {
            auto contraction = t1.contract(t2, ContractionPairs<IndexI, IndexJ>());
            if constexpr (IsIdentity(perm))
                return contraction;
            else
                return contraction.shuffle(perm);
        }
    }

#pragma endregion
    //
    //
}
