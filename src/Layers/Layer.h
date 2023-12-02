#pragma once

#include<pch.h>
#include"Events/Events.h"
#include"Utilities/Memory.h"

namespace SCESS
{
	class Framebuffer;
	class Camera;

	class Layer
	{
	public:
		Layer(const std::string_view name = "Layer");
		virtual ~Layer();

		virtual void OnAttach() {}
		virtual void OnDetach() {}
		virtual void OnUpdate(const sRef<SCESS::Framebuffer>& screenBuffer, const Camera& cam,
			const glm::mat4& projectionMatrix, const glm::vec2& nearFarPlanes) 
			{}
		virtual void OnImGuiRender() {}
		virtual void OnEvent(Event& e) {}

		void SetRenderStyle(GLenum renderMode) { mRenderStyle = renderMode; }
		GLenum GetRenderStyle() const { return mRenderStyle; }

		virtual void Resize(unsigned int width, unsigned int height) {};

		inline const std::string& GetName() const { return m_Name; }
	protected:

		GLenum mRenderStyle = GL_FILL;
		std::string m_Name;
	};
}
