#pragma once

#include <algorithm>
#include "glm/vec2.hpp"
#include "glm/matrix.hpp"
#include "GraphicsEngine/GEng_API.h"

namespace GEng {
    GENG_API class Camera2D {
    public:
        Camera2D() = default;

        explicit Camera2D(const glm::vec2 nearFar): mProjectionMatrix(1.0f), mViewMatrix(1.0f), mNearFar(nearFar) {
        }

        Camera2D(const Camera2D &cam) : mViewMatrix(cam.mViewMatrix),
                                        mProjectionMatrix(cam.mProjectionMatrix),
                                        mPos(cam.mPos),
                                        mZoom(cam.mZoom),
                                        mNearFar(cam.mNearFar) {
        }

        Camera2D &operator=(const Camera2D &cam) {
            mViewMatrix = cam.mViewMatrix;
            mPos = cam.mPos;
            mZoom = cam.mZoom;

            return *this;
        }

        void UpdateMouse(float mouseDX,
                         float mouseDY,
                         int buttonClicked,
                         int modKeys = 0);

        void UpdateZoom(float mouseWheelD);

        void MoveTo(const glm::vec2 &newPos, const float newZoom) {
            mPos = newPos;
            mZoom = newZoom;
        }

        void UpdateCamera(float deltaTime);

        const glm::mat4 &GetViewMatrix() const {
            return mViewMatrix;
        }

        const glm::mat4 &GetProjectionMatrix() const { return mProjectionMatrix; }

        void SetPos(const glm::vec2 &pos) {
            mVel = {0.0f, 0.0f};
            mPos = pos;
        }

        const glm::vec2 &GetPos() const { return mPos; }

        const glm::vec2 &GetVelocity() const { return mVel; }
        const float &GetZoomVelocity() const { return mZoomVel; }

        const void SetZoom(float z) {
            mZoomVel = 0.0f;
            mZoom = std::clamp(z, 0.05f, 100.0f);
        }

        const float &GetZoom() const { return mZoom; }

        const glm::vec2 &GetNearFarPlanes() const { return mNearFar; }
        void SetNearFarPlanes(const glm::vec2 nearFar) { mNearFar = nearFar; }

        void RecalculateViewMatrix();
        void RecalculateOrthographicProjection(const glm::vec2 screenDimensions);

    private:
        glm::mat4 mProjectionMatrix;
        glm::mat4 mViewMatrix;
        glm::vec2 mPos = {0.0f, 0.0f};
        glm::vec2 mVel = {0, 0};
        glm::vec2 mNearFar = {0.01f, 1000.0f};

        float mZoom = 1.0f;
        float mZoomVel = 0.0f;
    };
}
