#include<pch.h>
//--
#include<filesystem>
//--
#include"ShaderProgram.h"
#include"Platform/API Core.h"
#include"Platform/OpenGL/OpenGLShaderProgram.h"
#include"Utilities/StringOperations.h"
#include"Utilities/IO.h"
#include"Utilities/Memory.h"
#include"Utilities/StringOperations.h"
//
//

namespace SCESS
{
	static bool removeContinue = false;

	void ShaderProgram::Preprocessor::PreprocessShaderType(
		const std::string combinedSourceCode, std::string* vertexOutput,
		std::string* fragmentOutput, std::string* geometryOutput)
	{
		std::string shaderHeader;
		if (size_t beginIndex = combinedSourceCode.find("#begin"); beginIndex != combinedSourceCode.npos)
		{
			shaderHeader = std::move(combinedSourceCode.substr(0, beginIndex));
		}

		auto EvaluateShader = [combinedSourceCode](const char* beginStr, const char* endStr, std::string& output) -> void
			{
				size_t beginIndex = combinedSourceCode.find(beginStr);
				size_t endIndex = combinedSourceCode.find(endStr);

				if (beginIndex != combinedSourceCode.npos)
				{
					endIndex = endIndex != combinedSourceCode.npos ? endIndex : combinedSourceCode.size() - 1;
					size_t beginNewLineIndex = combinedSourceCode.find('\n', beginIndex) + 1;
					output = combinedSourceCode.substr(beginNewLineIndex, endIndex - beginNewLineIndex);
				}
			};

		if (vertexOutput != nullptr)
		{
			EvaluateShader("#begin(vShader)", "#end(vShader)", *vertexOutput);
			(*vertexOutput) = shaderHeader + (*vertexOutput);
		}
		if (fragmentOutput != nullptr)
		{
			EvaluateShader("#begin(fShader)", "#end(fShader)", *fragmentOutput);
			(*fragmentOutput) = shaderHeader + (*fragmentOutput);
		}
		if (geometryOutput != nullptr)
		{
			EvaluateShader("#begin(gShader)", "#end(gShader)", *geometryOutput);
			(*geometryOutput) = shaderHeader + (*geometryOutput);
		}
	}

	void ShaderProgram::Preprocessor::Preprocess(const std::string_view shaderFilePath, std::string& sourceCode, ShaderKeyMap* insertionKeys)
	{
		std::vector<std::string> onceIncludeHeaders;

		std::vector<std::string> lines = StringOperations::Split(sourceCode, '\n');
		sourceCode.clear();

		for (std::string& currentLine : lines)
		{
			EvaluatePreprocessorDirectives(shaderFilePath, onceIncludeHeaders, currentLine, sourceCode, "", insertionKeys);
		}
		removeContinue = false;
		sourceCode.erase(sourceCode.size() - 1, 1);
		//sourceCode.SpliceAssign(sourceCode.GetLength() - 1, 1);

		//sourceCode.TrimToCapacity();
	}

	void ShaderProgram::Preprocessor::EvaluatePreprocessorDirectives(const std::string_view shaderFilePath, std::vector<std::string>& onceIncludeHeaders, std::string& currentLine, std::string& sourceCode, const std::string_view lastInclude,
		ShaderKeyMap* insertionKeys)
	{
		if (currentLine.find("#pragma endremove") != currentLine.npos)
		{
			removeContinue = false;
			return;
		}
		else if (currentLine.find("#pragma remove") != currentLine.npos)
		{
			removeContinue = true;
			return;
		}
		if (removeContinue) return;

		size_t pIndex = currentLine.find("#pragma");
		if (pIndex != currentLine.npos)
		{
			if (size_t insertIndex = currentLine.find("insert", pIndex); insertIndex != currentLine.npos)
			{
				size_t parBeginIndex = currentLine.find('(', insertIndex + 6);
				if (parBeginIndex == currentLine.npos)
				{
					UNI_WARN(SCESS::EvaluatePreprocessorDirectives, "Please begin #include directive file path with < or \"");
					return;
				}
				size_t parEndIndex = currentLine.find(')', insertIndex + 6);
				if (parEndIndex == currentLine.npos)
				{
					UNI_WARN(SCESS::ShaderProgram::Preprocessor::EvaluatePreprocessorDirectives, "Please begin #include directive file path with < or \"");
					return;
				}
				std::string insertKey = currentLine.substr(parBeginIndex + 1, parEndIndex - parBeginIndex - 1);
				if (insertionKeys == nullptr)
				{
					UNI_WARN(SCESS::ShaderProgram::Preprocessor::EvaluatePreprocessorDirectives, "Insertion keys map was nullptr, thus #pragma insert ");
					//sourceCode += currentLine.Substring(0, pIndex) + currentLine.Substring(parEndIndex, currentLine.GetLength() - parEndIndex) + newLineChar;
					currentLine.erase(pIndex, (parEndIndex + 1) - pIndex);
					sourceCode += currentLine + '\n';
					return;
				}

				auto insertKVPair = insertionKeys->find(insertKey);
				if (insertKVPair == insertionKeys->end())
				{
					UNI_WARN(SCESS::ShaderProgram::Preprocessor::EvaluatePreprocessorDirectives, ("#pragma insert(" + insertKey + ") could not be evaluated, as key map contains no member '" + insertKey + '\'').c_str());
					//sourceCode += currentLine.Substring(0, pIndex) + currentLine.Substring(parEndIndex, currentLine.GetLength() - parEndIndex) + newLineChar;
					currentLine.erase(pIndex, (parEndIndex + 1) - pIndex);
					sourceCode += currentLine + '\n';
					return;
				}
				sourceCode += currentLine.substr(0, pIndex) + insertKVPair->second + currentLine.substr(parEndIndex + 1, currentLine.size() - parEndIndex - 1) + '\n';
			}
		}

		size_t includeIndex = currentLine.find("#include");
		if (includeIndex != currentLine.npos)
		{
			size_t beginIndex = currentLine.find('<', includeIndex + 8);
			beginIndex = beginIndex == currentLine.npos ? currentLine.find('\"', includeIndex + 8) : beginIndex;
			if (beginIndex == currentLine.npos)
			{
				UNI_FATAL(SCESS::ShaderProgram::Preprocessor::EvaluatePreprocessorDirectives, "Please begin #include directive file path with < or \"");
				return;
			}

			size_t endIndex = currentLine.find('>', beginIndex + 1);
			endIndex = endIndex == currentLine.npos ? currentLine.find('\"', beginIndex + 1) : endIndex;
			if (endIndex == currentLine.npos)
			{
				UNI_FATAL(SCESS::ShaderProgram::Preprocessor::EvaluatePreprocessorDirectives, "Please end #include directive file path with > or \"");
				return;
			}

			size_t commentIndex = currentLine.find("//", 0);
			if (commentIndex < includeIndex || commentIndex < endIndex)
			{
				sourceCode += (currentLine + '\n');
				return;
			}


			std::string includePath = currentLine.substr(beginIndex + 1, endIndex - beginIndex - 1);
			StringOperations::Replace(includePath, "/", "\\");
			if (includePath.rfind("..\\", 0) == 0)
			{
				includePath.insert(0, "..\\");
			}
			std::string filePath = std::operator+(std::string(IO::Dir::GetDirectoryFromPath(shaderFilePath)), includePath);
			if (filePath == lastInclude)
			{
				UNI_FATAL(EvaluatePreprocessorDirectives, "Infinite inclusion loop. File " + filePath + " includes itself.");
				return;
			}
			if (filePath.size() == 0)
			{
				UNI_FATAL(SCESS::ShaderProgram::Preprocessor::EvaluatePreprocessorDirectives, "No file was given when evaluating #include directive");
				return;
			}


			if (std::find(onceIncludeHeaders.begin(), onceIncludeHeaders.end(), filePath) == onceIncludeHeaders.end())
			{
				std::string fileData = std::move(IO::ReadTextFile(filePath));
				pIndex = fileData.find("#pragma");
				if (pIndex != currentLine.npos)
				{
					size_t onceIndex = fileData.find("once", pIndex);
					if (onceIndex != currentLine.npos)
					{
						onceIncludeHeaders.push_back(filePath);
						fileData.erase(pIndex, (onceIndex + 4) - pIndex);
					}
				}
				if (fileData.find("#include") != fileData.npos || pIndex != currentLine.npos)
				{
					std::vector<std::string> fileDataLines = StringOperations::Split(fileData, '\n');
					fileData.clear();

					for (std::string& line : fileDataLines)
					{
						EvaluatePreprocessorDirectives(shaderFilePath, onceIncludeHeaders, line, fileData, filePath, insertionKeys);
					}
				}
				sourceCode += (currentLine.substr(0, includeIndex) + fileData + currentLine.substr(endIndex + 1, currentLine.size()) + '\n');
				return;
			}
			else
			{
				sourceCode += currentLine.substr(0, includeIndex) + currentLine.substr(endIndex + 1, currentLine.size()) + '\n';
				return;
			}
		}

		if (pIndex == currentLine.npos && includeIndex == currentLine.npos) sourceCode += currentLine + '\n';
	}

	sRef<ShaderProgram> ShaderProgram::Create(const std::string_view combinedSourceCode)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
		{
			sRef<OpenGLShaderProgram> program = make <OpenGLShaderProgram>();
			program->CompileCombinedSource(combinedSourceCode);
			return program;
		}
		}
	}

	sRef<ShaderProgram> ShaderProgram::Create(const std::string_view vertexSource, const std::string_view fragmentSource)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
		{
			sRef<OpenGLShaderProgram> program = make <OpenGLShaderProgram>();
			program->CompileSource(vertexSource, fragmentSource);
			return program;
		}
		}
	}

	sRef<ShaderProgram> ShaderProgram::CreateFromFile(const std::string_view combinedShaderFile)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
		{
			sRef<OpenGLShaderProgram> program = make <OpenGLShaderProgram>();
			program->CompileCombinedFile(combinedShaderFile);
			return program;
		}
		}
	}

	sRef<ShaderProgram> ShaderProgram::CreateFromFile(const std::string_view vertexFile, const std::string_view fragmentFile)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
		{
			sRef<OpenGLShaderProgram> program = make <OpenGLShaderProgram>();
			program->CompileFile(vertexFile, fragmentFile);
			return program;
		}
		}
	}
}
