#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Events/Events.h"
#include <glm/glm.hpp>

namespace GEng
{
	// Size vector of dimension 2
	using svec2 = glm::vec<2, size_t>;

	GENG_API enum class MouseButton
	{
		None = -1,
		Right = GLFW_MOUSE_BUTTON_RIGHT,
		Centre = GLFW_MOUSE_BUTTON_MIDDLE,
		Left = GLFW_MOUSE_BUTTON_LEFT
	};

	GENG_API class MouseButtonPressedEvent : public Event
	{
	public:
		MouseButtonPressedEvent(MouseButton buttonIn) : button(buttonIn) {}

		EVENT_CLASS_FUNCTIONS(MouseButtonPressed)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryInput |
								 EventCategory::EventCategoryMouse |
								 EventCategory::EventCategoryMouseButton)

		inline MouseButton GetMouseButton() const { return button; }

	protected:
		MouseButton button;
	};

	GENG_API class MouseButtonReleasedEvent : public Event
	{
	public:
		MouseButtonReleasedEvent(MouseButton buttonIn) : button(buttonIn) {}

		EVENT_CLASS_FUNCTIONS(MouseButtonReleased)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryInput |
								 EventCategory::EventCategoryMouse |
								 EventCategory::EventCategoryMouseButton)

		inline MouseButton GetMouseButton() const { return button; }

	protected:
		MouseButton button;
	};

	GENG_API class MouseMoveEvent : public Event
	{
	public:
		MouseMoveEvent(glm::dvec2 pos) : mMousePos(pos) {}

		EVENT_CLASS_FUNCTIONS(MouseMove)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryInput | EventCategory::EventCategoryMouse)

		// MouseButton GetMouseButtonPressed() const { return buttonPressed; }
		glm::dvec2 GetMousePos() const { return mMousePos; }

	protected:
		glm::dvec2 mMousePos;
		// MouseButton buttonPressed;
	};

	GENG_API class MouseScrollEvent : public Event
	{
	public:
		MouseScrollEvent(float xOffIn, float yOffIn) : xOffset(xOffIn), yOffset(yOffIn) {}

		EVENT_CLASS_FUNCTIONS(MouseWheelScroll)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryInput | EventCategory::EventCategoryMouse)

		inline float GetXOffset() const { return xOffset; }
		inline float GetYOffset() const { return yOffset; }

	protected:
		float xOffset = 0.0f;
		float yOffset = 0.0f;
	};

}