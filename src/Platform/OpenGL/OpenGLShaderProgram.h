#pragma once
#include<pch.h>
#include<string_view>
#include<glm/glm.hpp>
#include"3D/ShaderProgram.h"

namespace SCESS
{
	class OpenGLShaderProgram : public ShaderProgram
	{
	public:
		using UniformLocMap = std::unordered_map<std::string_view, int>;

		OpenGLShaderProgram();

		void CompileSource(const std::string_view vertexSourceIn, const std::string_view fragmentSourceIn,
			const std::string_view vFilePath = "", const std::string_view fFilePath = "");
		void CompileCombinedSource(const std::string_view combinedSource, const std::string_view filePath = "");
		void CompileFile(const std::string_view vFile, const std::string_view fFile);
		void CompileCombinedFile(const std::string_view combinedShaderFile);

		bool Link();

		virtual ~OpenGLShaderProgram();

		virtual void Bind() override;
		virtual void Unbind() override;

		virtual void UploadF   (const std::string_view uniformName, const float val)      override;
		virtual void UploadVecF(const std::string_view uniformName, const glm::vec2& val) override;
		virtual void UploadVecF(const std::string_view uniformName, const glm::vec3& val) override;
		virtual void UploadVecF(const std::string_view uniformName, const glm::vec4& val) override;
		virtual void UploadMatF(const std::string_view uniformName, const glm::mat4& val) override;
		virtual void UploadMatF(const std::string_view uniformName, const glm::mat3& val) override;

		virtual void UploadI   (const std::string_view uniformName, const int val)         override;
		virtual void UploadVecI(const std::string_view uniformName, const glm::ivec2& val) override;
		virtual void UploadVecI(const std::string_view uniformName, const glm::ivec3& val) override;
		virtual void UploadVecI(const std::string_view uniformName, const glm::ivec4& val) override;

	protected:
		int QueryLocation(const std::string_view uniformName);


	private:
		UniformLocMap m_UniformMap;
		GLid m_ProgramID = 0;

#ifdef CONFIG_DEBUG
		std::string mFileName;
#endif // CONFIG_DEBUG
	};
}