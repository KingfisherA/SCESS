#include<pch.h>
#include "OpenGLShaderProgram.h"
#include"Utilities/Assert.h"
#include"Utilities/StringOperations.h"
#include"Utilities/IO.h"
//
//


namespace SCESS
{
	OpenGLShaderProgram::OpenGLShaderProgram()
	{
		m_ProgramID = glCreateProgram();
	}

	void OpenGLShaderProgram::CompileSource(const std::string_view vertexSourceIn, const std::string_view fragmentSourceIn,
		const std::string_view vFilePath, const std::string_view fFilePath)
	{
		std::string vertexSource = vertexSourceIn.data();
		std::string fragmentSource = fragmentSourceIn.data();
		Preprocessor::Preprocess(vFilePath, vertexSource, nullptr);
		Preprocessor::Preprocess(fFilePath, fragmentSource, nullptr);

		
		auto CheckShaderCompilation = [this](GLuint shaderID) {
			GLint compileError = 0;
			glGetShaderiv(shaderID, GL_COMPILE_STATUS, &compileError);
			if (compileError == GL_FALSE)
			{
				GLint logLength = 0;
				glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &logLength);

				std::string logString(logLength, ' ');
				glGetShaderInfoLog(shaderID, logLength, &logLength, &logString[0]);
#ifdef CONFIG_DEBUG
				
				logString = std::format("Shader compilation failed.\n Shader file path: {} \n Open GL log\n {}", mFileName, logString);
#else
				logString = "Shader compilation failed: " + logString;
#endif // CONFIG_DEBUG


				UNI_FATAL(OpenGLShaderProgram::CompileSource, logString);

				return false;
			}
			return true;
		};

		const GLchar* vertexSourcePtr = vertexSource.data();
		GLuint vertexShaderID = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vertexShaderID, 1, &vertexSourcePtr, 0);
		glCompileShader(vertexShaderID);
		if (!CheckShaderCompilation(vertexShaderID))
		{
			glDeleteShader(vertexShaderID);
			return;
		}

		const GLchar* fragmentSourcePtr = fragmentSource.data();
		GLuint fragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fragmentShaderID, 1, &fragmentSourcePtr, 0);
		glCompileShader(fragmentShaderID);
		if (!CheckShaderCompilation(fragmentShaderID))
		{
			glDeleteShader(vertexShaderID);
			glDeleteShader(fragmentShaderID);
			return;
		}

		glAttachShader(m_ProgramID, vertexShaderID);
		glAttachShader(m_ProgramID, fragmentShaderID);

		Link();

		glDetachShader(m_ProgramID, vertexShaderID);
		glDetachShader(m_ProgramID, fragmentShaderID);
		glDeleteShader(vertexShaderID);
		glDeleteShader(fragmentShaderID);
	}

	void OpenGLShaderProgram::CompileCombinedSource(const std::string_view combinedSource,
		const std::string_view filePath)
	{
		//std::string vSource = combinedSource.data();
		//size_t location = vSource.find("#version");
		//location = location == vSource.npos ? 0 : vSource.find('\n', location) + 1;
		//vSource.insert(location, "#define VERTEX_SHADER\n", 22);

		//std::string fSource = combinedSource.data();
		//location = fSource.find("#version");
		//location = location == fSource.npos ? 0 : fSource.find('\n', location) + 1;
		//fSource.insert(location, "#define FRAGMENT_SHADER\n", 24);

		std::string vSource, fSource;
		Preprocessor::PreprocessShaderType(combinedSource.data(), &vSource, &fSource, nullptr);

		this->CompileSource(vSource, fSource, filePath, filePath);
	}

	void OpenGLShaderProgram::CompileFile(const std::string_view vFile, const std::string_view fFile)
	{
		std::string vSourceCode = std::move(IO::ReadTextFile(vFile));
		std::string fSourceCode = std::move(IO::ReadTextFile(fFile));
		
#ifdef CONFIG_DEBUG
		mFileName = vFile;
#endif // CONFIG_DEBUG


		CompileSource(vSourceCode, fSourceCode, vFile, fFile);
	}

	void OpenGLShaderProgram::CompileCombinedFile(const std::string_view combinedShaderFile)
	{
		std::string sourceCode = std::move(IO::ReadTextFile(combinedShaderFile));

#ifdef CONFIG_DEBUG
		mFileName = combinedShaderFile;
#endif // CONFIG_DEBUG

		CompileCombinedSource(sourceCode, combinedShaderFile);
	}

	bool OpenGLShaderProgram::Link()
	{
		glLinkProgram(m_ProgramID);

		GLint linkStatus = 0;
		glGetProgramiv(m_ProgramID, GL_LINK_STATUS, &linkStatus);
		if (linkStatus == GL_FALSE)
		{
			GLint infoLogLength = 0;
			glGetProgramiv(m_ProgramID, GL_INFO_LOG_LENGTH, &infoLogLength);

			std::string infoLog(infoLogLength, ' ');
			glGetProgramInfoLog(m_ProgramID, infoLogLength, &infoLogLength, &infoLog[0]);
#ifdef CONFIG_DEBUG
			infoLog = std::format("Program Linking Failed.\n Shader file path: {} \n Open GL log\n {}", mFileName, infoLog);
#else
			infoLog = "Shader compilation failed: " + infoLog;
#endif // CONFIG_DEBUG


			UNI_FATAL(OpenGLShaderProgram::Link, infoLog);
		}

		return true;
	}

	OpenGLShaderProgram::~OpenGLShaderProgram()
	{
		glDeleteProgram(m_ProgramID);
	}

	void OpenGLShaderProgram::Bind()
	{
		glUseProgram(m_ProgramID);
	}

	void OpenGLShaderProgram::Unbind()
	{
		glUseProgram(0);
	}

	void OpenGLShaderProgram::UploadF   (const std::string_view uniformName, const float val)
	{
		glUniform1f(QueryLocation(uniformName), val);
	}

	void OpenGLShaderProgram::UploadVecF(const std::string_view uniformName, const glm::vec2& val)
	{
		glUniform2f(QueryLocation(uniformName), val.x, val.y);
	}

	void OpenGLShaderProgram::UploadVecF(const std::string_view uniformName, const glm::vec3& val)
	{
		glUniform3f(QueryLocation(uniformName), val.x, val.y, val.z);
	}

	void OpenGLShaderProgram::UploadVecF(const std::string_view uniformName, const glm::vec4& val)
	{
		glUniform4f(QueryLocation(uniformName), val.x, val.y, val.z, val.w);
	}

	void OpenGLShaderProgram::UploadMatF(const std::string_view uniformName, const glm::mat4& val)
	{
		glUniformMatrix4fv(QueryLocation(uniformName), 1, GL_FALSE, glm::value_ptr(val));
	}

	void OpenGLShaderProgram::UploadMatF(const std::string_view uniformName, const glm::mat3& val)
	{
		glUniformMatrix3fv(QueryLocation(uniformName), 1, GL_FALSE, glm::value_ptr(val));
	}

	void OpenGLShaderProgram::UploadI   (const std::string_view uniformName, const int val)
	{
		glUniform1i(QueryLocation(uniformName), val);
	}

	void OpenGLShaderProgram::UploadVecI(const std::string_view uniformName, const glm::ivec2& val)
	{
		glUniform2i(QueryLocation(uniformName), val.x, val.y);
	}

	void OpenGLShaderProgram::UploadVecI(const std::string_view uniformName, const glm::ivec3& val)
	{
		glUniform3i(QueryLocation(uniformName), val.x, val.y, val.z);
	}

	void OpenGLShaderProgram::UploadVecI(const std::string_view uniformName, const glm::ivec4& val)
	{
		glUniform4i(QueryLocation(uniformName), val.x, val.y, val.z, val.w);
	}

	int OpenGLShaderProgram::QueryLocation(const std::string_view uniformName)
	{
		auto it = m_UniformMap.find(uniformName);
		if (it != m_UniformMap.end())
		{
			return it->second;
		}

		int location = glGetUniformLocation(m_ProgramID, uniformName.data());
		if (location < 0)
		{
			std::cerr << "Uniform location for variable \"" << uniformName << "\" could not be found \n";
		}

		m_UniformMap.emplace(uniformName, location);

		return location;
	}
}