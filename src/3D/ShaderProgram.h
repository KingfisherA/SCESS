#pragma once
#include<string_view>

#include"Utilities/Memory.h"

namespace SCESS
{
	class ShaderProgram
	{
	public:
		using ShaderKeyMap = std::unordered_map<std::string, std::string>;

		virtual ~ShaderProgram() {}

		virtual void Bind() = 0;
		virtual void Unbind() = 0;

		virtual void UploadF   (const std::string_view name, const float val)      =0;
		virtual void UploadVecF(const std::string_view name, const glm::vec2& val) =0;
		virtual void UploadVecF(const std::string_view name, const glm::vec3& val) =0;
		virtual void UploadVecF(const std::string_view name, const glm::vec4& val) =0;
		virtual void UploadMatF(const std::string_view name, const glm::mat4& val) =0;
		virtual void UploadMatF(const std::string_view name, const glm::mat3& val) =0;

		virtual void UploadI   (const std::string_view name, const int val)         =0;
		virtual void UploadVecI(const std::string_view name, const glm::ivec2& val) =0;
		virtual void UploadVecI(const std::string_view name, const glm::ivec3& val) =0;
		virtual void UploadVecI(const std::string_view name, const glm::ivec4& val) =0;

		static sRef<ShaderProgram> Create(const std::string_view combinedSourceCode);
		static sRef<ShaderProgram> Create(const std::string_view vertexSource, const std::string_view fragmentSource);
		static sRef<ShaderProgram> CreateFromFile(const std::string_view combinedShaderFile);
		static sRef<ShaderProgram> CreateFromFile(const std::string_view vertexFile, const std::string_view fragmentFile);

		struct Preprocessor
		{
		public:
			static void PreprocessShaderType(const std::string combinedSourceCode, std::string* vertexOutput,
				std::string* fragmentOutput, std::string* geometryOutput);
			static void Preprocess(const std::string_view shaderFilePath, std::string& sourceCode, ShaderKeyMap* insertionKeys);
		private:
			static void EvaluatePreprocessorDirectives(const std::string_view shaderFilePath, std::vector<std::string>& onceIncludeHeaders, std::string& currentLine, std::string& sourceCode, const std::string_view lastInclude,
				ShaderKeyMap* insertionKeys);
		};
};
}