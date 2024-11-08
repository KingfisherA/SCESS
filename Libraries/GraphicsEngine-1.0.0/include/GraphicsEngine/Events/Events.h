#pragma once

#include "GraphicsEngine/GEng_API.h"
#include <string>
#include <string_view>
#include <functional>

namespace GEng
{
	 class Event;
	 class EventDispatcher;

	GENG_API enum class EventType
	{
		None = 0,
		WindowClose,
		WindowResize,
		WindowFocus,
		WindowMoved,
		AppTick,
		AppUpdate,
		AppRender,
		KeyPressed,
		KeyReleased,
		KeyTyped,
		MouseButtonPressed,
		MouseButtonReleased,
		MouseMove,
		MouseWheelScroll
	};

#define BIT(nr) (1u << (nr))
	enum EventCategory
	{
		None = 0,
		EventCategoryApplication = BIT(0),
		EventCategoryInput = BIT(1),
		EventCategoryKeyboard = BIT(2),
		EventCategoryMouse = BIT(3),
		EventCategoryMouseButton = BIT(4),
		EventCategoryWindow = BIT(5)
	};

#define EVENT_CLASS_FUNCTIONS(type)                                             \
	 static EventType GetStaticType() { return EventType ::type; }               \
	 virtual EventType GetEventType() const override { return GetStaticType(); } \
	 virtual std::string_view GetName() const override { return std::string_view(#type); }

#define EVENT_CLASS_CATEGORY_FNS(cat) \
	virtual int GetCategoryFlags() const override { return cat; }

	GENG_API class Event
	{
		friend GENG_API class EventDispatcher;

	public:
		virtual EventType GetEventType() const = 0;
		virtual std::string_view GetName() const = 0;
		virtual int GetCategoryFlags() const = 0;
		virtual std::string ToString() const { return std::string(GetName()); }
		inline bool IsCategory(EventCategory cat) { return GetCategoryFlags() & cat; }
		inline void Handle() { m_Handled = true; }
		inline bool IsHandled() const { return m_Handled; }

	protected:
		bool m_Handled = false;
	};

	GENG_API class EventDispatcher
	{
	private:
		Event &m_Event;

	public:
		template <typename T>
		using EventFn = std::function<bool(T &)>;

		EventDispatcher(Event &event) : m_Event(event)
		{
		}

		template <typename T>
		bool Dispatch(EventFn<T> func)
		{
			if (m_Event.GetEventType() == T::GetStaticType() && !m_Event.m_Handled)
			{
				m_Event.m_Handled = func((T &)m_Event);
				return true;
			}
			return false;
		}

		const Event &GetEvent() const { return m_Event; }
	};

}
