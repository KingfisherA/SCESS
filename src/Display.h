#pragma once
#include"Events/Events.h"


namespace SCESS
{
	struct DisplayProps
	{
		std::string title;
		size_t width = 0, height = 0;
	};

	class Display
	{
	public:
		using EventCallbackFnc = std::function<void(Event&)>;

		virtual ~Display() {}

		virtual void Refresh() = 0;
		virtual size_t GetWidth() const = 0;
		virtual size_t GetHeight() const = 0;

		virtual void Show() const = 0;
		virtual void Hide() const = 0;
		virtual void Focus() const = 0;

		virtual bool CloseRequested() const = 0;
		virtual double GetTime() const = 0;

		virtual void BindEventCallback(const EventCallbackFnc& callback) = 0;

		virtual void* GetAPIDisplay() const = 0;

		static Display* Create(const DisplayProps& properties);
	};
}
