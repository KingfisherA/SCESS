#pragma once
#include"Layers/Layer.h"
#include"Events/Events.h"

namespace SCESS
{
	class ImGuiLayer : public Layer
	{
	public:
		ImGuiLayer();
		~ImGuiLayer();

		virtual void OnAttach() override;
		virtual void OnDetach() override;
		virtual void OnImGuiRender() override;

		virtual void OnEvent(Event& e) override;

		void RenderDockspace();

		void Begin();
		void End();

	private:
		void SetDarkThemeColors();
	};
}