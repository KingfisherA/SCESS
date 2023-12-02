#pragma once
#include"Events.h"
#include<glm/glm.hpp>

namespace SCESS
{
	//Size vector of dimension 2
	using svec2 = glm::vec<2, size_t>;

	enum class MouseButton
	{
		None = -1, 
		Right = GLFW_MOUSE_BUTTON_RIGHT, Centre = GLFW_MOUSE_BUTTON_MIDDLE, Left = GLFW_MOUSE_BUTTON_LEFT
	};

	class MouseButtonPressedEvent : public Event
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

	class MouseButtonReleasedEvent : public Event
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

	class MouseMoveEvent : public Event
	{
	public:
		MouseMoveEvent(glm::vec2 pos) : mousePos(pos) {}

		EVENT_CLASS_FUNCTIONS(MouseMove)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryInput | EventCategory::EventCategoryMouse)

		//MouseButton GetMouseButtonPressed() const { return buttonPressed; }
		svec2 GetMousePos() const { return mousePos; }

	protected:
		glm::vec2 mousePos = { 0,0 };
		//MouseButton buttonPressed;
	};

	class MouseScrollEvent : public Event
	{
	public:
		MouseScrollEvent(float xOffIn, float yOffIn) : xOffset(xOffIn), yOffset(yOffIn) {}

		EVENT_CLASS_FUNCTIONS(MouseWheelScroll)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryInput|EventCategory::EventCategoryMouse)

		inline float GetXOffset() const { return xOffset; }
		inline float GetYOffset() const { return yOffset; }

	protected:
		float xOffset = 0.0f;
		float yOffset = 0.0f;
	};

}