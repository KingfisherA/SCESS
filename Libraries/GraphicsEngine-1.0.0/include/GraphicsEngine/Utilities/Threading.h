#pragma once

#include "GraphicsEngine/GEng_API.h"
#include <thread>

namespace GEng {
    GENG_API class Thread {
    public:
        Thread() : mThread(), mThreadFinished(false), mThreadStarted(false) {
        }

        Thread(const Thread &) = delete;

        Thread(Thread &&t) noexcept : mThread(std::move(t.mThread)), mThreadFinished(t.mThreadFinished.load()),
                                      mThreadStarted(true) {
        }

        virtual ~Thread() {
            if (mThread.joinable()) {
                mThread.request_stop();
                mThread.join();
            }
        }

        Thread &operator=(Thread &&t) noexcept {
            Stop();
            mThread = std::move(t.mThread);
            mThreadFinished.store(t.mThreadFinished.load());
            mThreadStarted = t.mThreadStarted;

            return *this;
        }

        template<class Fn, class... Args>
        void Start(Fn &&Fx, Args &&... Ax) {
            Stop();
            mThreadFinished.store(false);
            mThreadStarted = true;
            mThread = std::move(std::jthread(std::forward<Fn>(Fx), std::forward<Args>(Ax)...));
        }

        void Stop() {
            if (!mThreadFinished.load() && mThread.joinable()) {
                mThreadFinished.store(true);
                mThread.request_stop();
                mThread.join();
            }
        }

        void Join() {
            if (!mThreadFinished.load() && mThread.joinable())
                mThread.join();
        }

        void Detach() {
            if (!mThreadFinished.load()) {
                mThread.detach();
            }
        }

        void MarkAsFinished() {
            mThreadFinished.store(true);
        }

        [[nodiscard]] bool ShouldStop() const { return mThread.get_stop_token().stop_requested(); }
        [[nodiscard]] bool IsRunning() const { return (!mThreadFinished.load()) && mThreadStarted; }
        [[nodiscard]] bool HasFinished() const { return mThreadFinished.load() && mThreadStarted; }

        [[nodiscard]] std::jthread &GetThread() { return mThread; }

    protected:
        std::jthread mThread;
        bool mThreadStarted = false;
        std::atomic<bool> mThreadFinished;
    };
}
