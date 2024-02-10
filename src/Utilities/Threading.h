#pragma once

#include<thread>

namespace SCESS
{
	class Thread
	{
	public:
		Thread() :mThread(nullptr), mThreadFinished(false) {  }
		Thread(const Thread&) = delete;
		Thread(Thread&& t) noexcept :mThread(std::move(t.mThread)) { mThreadFinished.store(t.mThreadFinished); };

		virtual ~Thread()
		{
			mThread->request_stop();
		}

		void operator=(Thread&& t) noexcept
		{
			Stop();
			mThread = std::move(t.mThread);
			mThreadFinished.store(t.mThreadFinished.load());
		}

		template <class _Fn, class... _Args>
		void Start(_Fn&& _Fx, _Args&&... _Ax)
		{
			Stop();
			mThreadFinished.store(false);
			mThread = std::make_unique<std::jthread>(std::forward<_Fn>(_Fx), std::forward<_Args>(_Ax)...);
		}

		void Stop()
		{
			if (mThread == nullptr)
				return;

			if (!mThreadFinished) {
				mThreadFinished.store(true);
				mThread->request_stop();
				mThread->join();
			}
		}

		void Join()
		{
			Stop();
		}

		void Detach()
		{
			if (mThread == nullptr)
				return;
			if (!mThreadFinished) {
				mThread->detach();
			}
		}

		void MarkAsFinished()
		{
			mThreadFinished.store(true);
		}

		const bool ShouldStop() const { return mThread->get_stop_token().stop_requested(); }
		const bool IsRunning() const { return (!mThreadFinished.load()) && mThread != nullptr; }
		const bool HasFinished() const { return mThreadFinished.load(); }

	protected:
		std::unique_ptr<std::jthread> mThread;
		std::atomic<bool> mThreadFinished;
	};
}