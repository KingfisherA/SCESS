#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Events/Events.h"

namespace GEng
{
	GENG_API class KeyEvent : public Event
	{
	public:
		KeyEvent(int keyCodeIn, int modsIn) : mKeyCode(keyCodeIn), mMods(modsIn) {}

		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryKeyboard | EventCategory::EventCategoryInput)

		const int GetKeyCode() const { return mKeyCode; }

		const int GetModifiers() const { return mMods; }

	protected:
		int mKeyCode, mMods;
	};

	GENG_API class KeyPressedEvent : public KeyEvent
	{
	public:
		KeyPressedEvent(int keyCode, int mods, unsigned int repeatCountIn) : KeyEvent(keyCode, mods), repeatCount(repeatCountIn) {}

		EVENT_CLASS_FUNCTIONS(KeyPressed)

		const unsigned int GetRepeatCount() const { return repeatCount; }

	protected:
		unsigned int repeatCount;
	};

	GENG_API class KeyReleasedEvent : public KeyEvent
	{
	public:
		KeyReleasedEvent(int keyCode, int mods) : KeyEvent(keyCode, mods) {}

		EVENT_CLASS_FUNCTIONS(KeyReleased)
	};

	GENG_API class KeyTypedEvent : public Event
	{
	public:
		KeyTypedEvent(unsigned int charCodeIn) : charCode(charCodeIn) {}

		EVENT_CLASS_FUNCTIONS(KeyTyped)
		EVENT_CLASS_CATEGORY_FNS(EventCategory::EventCategoryKeyboard | EventCategory::EventCategoryInput)

		const unsigned int GetCharCode() const { return charCode; }

	protected:
		unsigned int charCode;
	};
}