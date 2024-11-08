#pragma once
#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Events/Events.h"

namespace GEng
{

	GENG_API class AppRender : public Event
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