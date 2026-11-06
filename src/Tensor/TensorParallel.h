#pragma once

#ifndef EIGEN_USE_THREADS
#error "EIGEN_USE_THREADS must be defined (in every translation unit) to evaluate tensor expressions in parallel"
#endif

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <vector>
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
    //Calls fn(item, threadIndex) for every item in [0, nItems) on every hardware thread, handing the items out to
    //the threads in order. threadIndex is in [0, NumComputeThreads()) so it can index per-thread buffers. The first
    //exception thrown by fn stops the loop and is rethrown.
    template<class Fn>
    void ParallelForDynamic(size_t nItems, Fn &&fn) {
        const int nThreads = (int) std::min<size_t>(NumComputeThreads(), std::max<size_t>(nItems, 1));
        std::atomic<size_t> next{0};
        std::exception_ptr error;
        std::mutex errorMutex;

        auto worker = [&](int thread) {
            try {
                for (size_t item; (item = next.fetch_add(1)) < nItems;)
                    fn(item, thread);
            } catch (...) {
                std::lock_guard lock(errorMutex);
                if (!error)
                    error = std::current_exception();
                next = nItems;
            }
        };

        std::vector<std::thread> threads;
        threads.reserve(nThreads - 1);
        for (int t = 1; t < nThreads; t++)
            threads.emplace_back(worker, t);
        worker(0);
        for (auto &thread: threads)
            thread.join();

        if (error)
            std::rethrow_exception(error);
    }
    //
    //
}
