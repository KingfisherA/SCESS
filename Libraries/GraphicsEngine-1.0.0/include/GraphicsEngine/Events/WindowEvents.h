#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Events/Events.h"
#include <sstream>
#include <glm/glm.hpp>

namespace GEng
{
	// Size vector of dimension 2
	using svec2 = glm::vec<2, size_t>;

	GENG_API class WindowResizeEvent : public Event
	{
	public:
		WindowResizeEvent(size_t widthIn, size_t heightIn) : width(widthIn),
															 height(heightIn) {}

		EVENT_CLASS_FUNCTIONS(WindowResize)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryWindow)

		inline size_t GetWidth() const { return width; }
		inline size_t GetHeight() const { return height; }
		inline svec2 GetDimensions() const { return {width, height}; }

		virtual std::string ToString() const override
		{
			std::stringstream ss;
			ss << "Window Resize Event: New Width: " << width << " New Height: " << height << ' ';
			return ss.str();
		}

	protected:
		size_t width = 0, height = 0;
	};

	GENG_API class WindowMovedEvent : public Event
	{
	public:
		WindowMovedEvent(svec2 newPos) : winPos(newPos) {}

		EVENT_CLASS_FUNCTIONS(WindowMoved)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryWindow)

		inline svec2 GetPosition() const { return winPos; }

	protected:
		svec2 winPos = {0, 0};
	};

	GENG_API class WindowFocusEvent : public Event
	{
	public:
		WindowFocusEvent(bool lostFocusIn) : lostFocus(lostFocusIn) {}

		EVENT_CLASS_FUNCTIONS(WindowFocus)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryWindow)

		inline bool DidLoseFocus() const { return lostFocus; }

	protected:
		bool lostFocus = false;
	};

	GENG_API class WindowCloseEvent : public Event
	{
	public:
		WindowCloseEvent() {}
		EVENT_CLASS_FUNCTIONS(WindowClose)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryWindow)
	};
}
