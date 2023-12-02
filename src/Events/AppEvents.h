#pragma once
#include"Events.h"

namespace SCESS
{

	class AppRender : public Event
	{
	public:
		AppRender(size_t deltaTime) : m_deltaTime(deltaTime) {}

		EVENT_CLASS_FUNCTIONS(AppRender)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryApplication)

		size_t GetDeltaFrameTime() const { return m_deltaTime; }

	protected:
		size_t m_deltaTime;
	};
}