#pragma once

#ifndef EIGEN_USE_THREADS
#error "EIGEN_USE_THREADS must be defined (in every translation unit) to evaluate tensor expressions in parallel"
#endif

#include <algorithm>
#include <thread>
#include <type_traits>
#include <unsupported/Eigen/CXX11/Tensor>

//Multithreaded evaluation of Eigen tensor expressions (such as those returned by Einsum) on a process-wide pool:
//  Tensor4 t = Evaluate(Einsum<...>(a, b));      evaluate into a new tensor
//  Parallel(t) += Einsum<...>(a, b);             evaluate into an existing tensor of the right size
//Do not use these from inside other worker threads; the pool already uses every hardware thread.
namespace SCESS {
    //
    //
    inline int NumComputeThreads() {
        static const int nThreads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
        return nThreads;
    }
    //
    inline const Eigen::ThreadPoolDevice &ParallelDevice() {
        static Eigen::ThreadPool pool(NumComputeThreads());
        static const Eigen::ThreadPoolDevice device(&pool, NumComputeThreads());
        return device;
    }
    //
    //Assignments to the returned object (=, +=, -=) are evaluated on the thread pool. The tensor must already have
    //the dimensions of the expression.
    template<class TensorType>
    auto Parallel(TensorType &t) {
        return t.device(ParallelDevice());
    }
    //
    //Evaluates a tensor expression into a new tensor on the thread pool
    template<class Expr>
    auto Evaluate(const Expr &expr) {
        using Traits = Eigen::internal::traits<Expr>;
        using Scalar = std::remove_const_t<typename Traits::Scalar>;
        constexpr int rank = Traits::NumDimensions;
        constexpr int layout = Traits::Layout;

        const Eigen::TensorEvaluator<const Expr, Eigen::ThreadPoolDevice> evaluator(expr, ParallelDevice());
        Eigen::Tensor<Scalar, rank, layout> result(evaluator.dimensions());
        result.device(ParallelDevice()) = expr;
        return result;
    }
    //
    //
}
