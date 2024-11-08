#pragma once

#include "GraphicsEngine/GEng_API.h"
#include <glm/glm.hpp>

namespace GEng {
    GENG_API class Camera3D {
    public:
        explicit Camera3D(float pitchYawAcc = 10100.0f, float posAcc = 5000.0f, float zoomAcc = 5000.0f,
                          float decayFactor = 0.91f);
        //
        void UpdateMovement(const glm::vec2 mouseDPosNormalised, int buttonClicked, int key = 0, int modKeys = 0);
        void UpdateZoom(double mouseWheelD);
        //
        [[nodiscard]] float GetPitch() const { return mPitch; }
        [[nodiscard]] float GetYaw() const { return mYaw; }
        [[nodiscard]] float GetRoll() const { return mRoll; }
        //
        [[nodiscard]] glm::vec3 GetRotation() const { return {mPitch, mYaw, mRoll}; }
        //
        void UpdateCamera(float deltaTime);
        //
        [[nodiscard]] const glm::mat4 &GetViewMatrix() const { return mViewMatrix; }
        [[nodiscard]] const glm::mat4 &GetInverseViewMatrix() const { return mInvViewMatrix; }
        [[nodiscard]] const glm::vec3 &GetCentrePosition() const { return mCentrePos; }
        void SetCentrePosition(const glm::vec3 &pos) {
            mCentrePos = pos;
        }
        [[nodiscard]] const glm::vec3 &GetActualPosition() const { return mActualPos; }
        [[nodiscard]] const float &GetZoom() const { return mZoom; }
        void SetZoom(const float zoomIn) {
            mZoom = zoomIn;
        }
        //
        [[nodiscard]] const float &GetX() const { return mCentrePos.x; }
        [[nodiscard]] const float &GetY() const { return mCentrePos.y; }
        [[nodiscard]] const float &GetZ() const { return mCentrePos.z; }
        //
        [[nodiscard]] const glm::vec3 &GetVelocity() const { return mPosVel; }
        //
        [[nodiscard]] const float &GetPitchYawAcc() const {
            return mPitchYawAcc;
        }
        void SetPitchYawAcc(float v) {
            mPitchYawAcc = v;
        }
        [[nodiscard]] const float &GetPosAcc() const {
            return mPosAcc;
        }
        void SetPosAcc(float v) {
            mPosAcc = v;
        }
        [[nodiscard]] const float &GetZoomAcc() const {
            return mZoomAcc;
        }
        void SetZoomAcc(float v) {
            mZoomAcc = v;
        }
        [[nodiscard]] const float &GetDecayFactor() const {
            return mDecayFactor;
        }
        void SetDecayFactor(float v) {
            mDecayFactor = v;
        }
        //
    private:
        //
        void RecalculateActualPos();
        void RecalculateViewMatrix();
        //
        glm::mat4 mViewMatrix;
        glm::mat4 mInvViewMatrix{};
        //
        glm::vec3 mPosVel;
        //
        glm::vec3 mCentrePos;
        glm::vec3 mActualPos{};
        //
        float mPitchYawAcc, mPosAcc, mZoomAcc, mDecayFactor;
        //
        float mRoll = 0.0f; // Z axis
        //
        float mPitch = 0.0f; // X axis
        float mPitchVel = 0.0f;
        //
        float mYaw = 0.0f; // Y axis
        float mYawVel = 0.0f;
        //
        float mZoom = 5.0f;
        float mZoomVel = 0.0f;
        //
    };
}
