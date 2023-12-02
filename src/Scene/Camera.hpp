#pragma once
#include<pch.h>

namespace SCESS
{
	class Camera
	{
	public:
		Camera();
		//
		void UpdateMouse(float mouseDX, float mouseDY, int buttonClicked);
		//
		void UpdatePitch(float mouseDY, int buttonClicked);
		void UpdateYaw(float mouseDX, int buttonClicked);
		void UpdateZoom(double mouseWheelD);
		//
		void MoveTo(const glm::vec3& newPosition, const float targetZoomIn);
		void StopPath();
		//
		inline float GetPitch() const { return pitch; }
		inline float GetYaw() const { return yaw; }
		inline float GetRoll() const { return roll; }
		//
		inline glm::vec3 GetRotation() const { return glm::vec3(pitch, yaw, roll); }
		//
		void UpdateCamera(double deltaTime, int key);
		//
		inline const glm::mat4& GetViewMatrix() const { return viewMatrix; }
		inline const glm::mat4& GetInverseViewMatrix() const { return invViewMatrix; }
		inline const glm::mat4& GetRotationsMatrix() const { return rotationsMatrix; }
		inline const glm::vec3& GetCentrePosition() const { return centrePosition; }
		void SetCentrePosition(const glm::vec3& pos) { StopPath(); centrePosition = pos; }
		const glm::vec3& GetActualPosition() const { return actualPosition; }
		const float GetZoom() const { return zoom; }
		void SetZoom(const float zoomIn) { StopPath(); zoom = zoomIn; }
		//
		inline float GetX() const { return centrePosition.x; }
		inline float GetY() const { return centrePosition.y; }
		inline float GetZ() const { return centrePosition.z; }
		//
		inline const glm::vec3& GetVelocity() const { return posVelocity; }
		//
	private:
		//
		void UpdateKeyboard(int key);
		void RecalculateActualPos();
		void RecalculateViewMatrix();
		//
		glm::mat4 rotationsMatrix;
		glm::mat4 viewMatrix;
		glm::mat4 invViewMatrix;
		//
		bool isFollowingPath = false;
		float posAcceleration = 0.0f;
		glm::vec3 posVelocity;
		//
		glm::vec3 centrePosition;
		glm::vec3 targetPosition;
		glm::vec3 actualPosition;
		//
		float roll = 0.0f; //Z axis
		//
		float pitch = 0.0f; //X axis
		float pitchVelocity = 0.0f;
		float pitchAcceleration = 0.0f;
		//
		float yaw = 0.0f; //Y axis
		float yawAcceleration = 0.0f;
		float yawVelocity = 0.0f;
		//
		float zoom = 5.0f;
		float zoomVelocity = 0.0f;
		float zoomAcceleration = 0.0f;
		float targetZoom = 0.0f;
	};
}