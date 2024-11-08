#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Events/Events.h"
#include "GraphicsEngine/Utilities/Memory.h"
#include "GraphicsEngine/Platform/OpenGL/OpenGLContext.h"

#include <glm/matrix.hpp>

namespace GEng
{
	class Framebuffer;
	class Camera3D;

	GENG_API class Layer
	{
	public:
		Layer(const std::string_view name = "Layer");
		virtual ~Layer();

		virtual void OnInit() {}
		virtual void OnDestroy() {}
		virtual void OnUpdate(const sRef<GEng::Framebuffer> &screenBuffer, const Camera3D &cam,
							  const glm::mat4 &projectionMatrix, const glm::vec2 &nearFarPlanes)
		{
		}
		virtual void OnImGuiRender() {}
		virtual void OnEvent(Event &e) {}

		void SetRenderStyle(const GPUContext::PolygonMode renderMode) { mPolygonMode = renderMode; }
		[[nodiscard]] GPUContext::PolygonMode GetRenderStyle() const { return mPolygonMode; }

		virtual void Resize(unsigned int width, unsigned int height) {};

		inline const std::string &GetName() const { return m_Name; }

	protected:
		GPUContext::PolygonMode mPolygonMode = GPUContext::PolygonMode::Solid;
		std::string m_Name;
	};
}
