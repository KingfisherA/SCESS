#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Layers/Layer.h"

namespace GEng
{
	GENG_API class LayerStack
	{
	public:
		using LayerVec = std::vector<Layer *>;

		LayerStack();
		~LayerStack();

		void PushLayer(Layer *layer);
		void PushOverlay(Layer *overlay);
		void PopLayer(Layer *layer);
		void PopOverlay(Layer *overlay);

		LayerVec::reverse_iterator rbegin() { return m_Layers.rbegin(); }
		LayerVec::reverse_iterator rend() { return m_Layers.rend(); }
		LayerVec::iterator begin() { return m_Layers.begin(); }
		LayerVec::iterator end() { return m_Layers.end(); }

	private:
		LayerVec m_Layers;
		size_t m_LayerInsertIndex;
	};
}
