#pragma once
//Adapted from https://github.com/nihaljn/marching-cubes
#include<glm/vec3.hpp>
#include<utility>
#include<condition_variable>
#include<mutex>
#include<memory>
#include<vector>
#include<functional>

#include "Dahlia/Tools/Threading.h"

namespace SCESS::Viewer {
    class ScalarField;
}

//
//
namespace SCESS::Viewer {
    //
    //
    struct GridCell {
        glm::vec3 vertex[8];
        float value[8];
    };
    //
    template<typename T>
    struct Triangle {
        typedef glm::vec<3, T, glm::packed_highp> VecT;

        Triangle() = default;

        Triangle(VecT vertex1, VecT vertex2, VecT vertex3) {
            mVertices[0] = vertex1;
            mVertices[1] = vertex2;
            mVertices[2] = vertex3;
        }

        Triangle(const Triangle &tri) {
            mVertices[0] = tri.mVertices[0];
            mVertices[1] = tri.mVertices[1];
            mVertices[2] = tri.mVertices[2];
        }

        Triangle &operator=(const Triangle &tri) {
            if (this == &tri)
                return *this;

            mVertices[0] = tri.mVertices[0];
            mVertices[1] = tri.mVertices[1];
            mVertices[2] = tri.mVertices[2];
            return *this;
        }

        VecT &operator[](unsigned int index) {
            return mVertices[index];
        }

        const VecT &operator[](unsigned int index) const {
            return mVertices[index];
        }

    public:
        VecT mVertices[3];
    };
    //
    //
    //
    class MarchingCubes {
    public:
        /// edgeToVertices[i] = {a, b} => edge i joins vertices a and b
        static const std::pair<int, int> edgeToVertices[12];

        /// edgeTable[i] is a 12 bit number; i is a cubeIndex
        /// edgeTable[i][j] = 1 if isosurface intersects edge j for cubeIndex i
        static const int edgeTable[256];

        /// triangleTable[i] is a list of edges forming triangles for cubeIndex i
        static const int triangleTable[256][16];

        /**
         * Given a `cell`, calculate its cube index
         * The cube index is an 8-bit encoding. Each bit represents a vertex. `index[i]` is the ith bit
         * If the value at the ith vertex is < isovalue, `index[i]` = 1. Else, `index[i]` = 0
         */
        int CalculateCubeIndex(GridCell &cell, float isovalue);

        /// Find the glm::vec3 between `v1` and `v2` where the functional value = `isovalue`
        glm::vec3 Interpolate(glm::vec3 &v1, float val1, glm::vec3 &v2, float val2, float isovalue);

        /// Returns all intersection coordinates of a cell with the isosurface
        /// (Calls `interpolate()`)
        void GenerateIntersectionCoordinates(GridCell &cell, float isovalue);

        [[nodiscard]] const glm::vec3 *GetIntersectionCoordinates() const { return mIntersections; }

        /// Given `cubeIndex`, get the edge table entry and using `intersections`, make all triangles
        void TriangulateIntersections(std::vector<Triangle<float> > &triangles, int cubeIndex);

        /// Get triangles of a single cell
        void TriangulateCell(std::vector<Triangle<float> > &triangles,
                             GridCell &cell, float isovalue);

        /// Triangulate a scalar field represented by `scalarFunction`. `isovalue` should be used for isovalue computation
        void TriangulateField(
            const ScalarField &field,
            float isovalue, std::vector<Triangle<float> > &trianglesOut, const DH::Thread *thread = nullptr);

    private:
        glm::vec3 mIntersections[12] = {};
    };
    //
    //
    //
    class ScalarField {
    public:
        typedef std::vector<std::vector<std::vector<float> > > TensorTy;

        ScalarField() = default;

        explicit ScalarField(const TensorTy &t) : mField(t) {
        }

        ScalarField(const ScalarField &s) : mField(s.mField) {
        }

        ScalarField(ScalarField &&s) noexcept : mField(std::move(s.mField)) {
        }

        void Initilise(size_t dimensions,
                       const std::function<float(float, float, float)> &func,
                       float isosurfaceVol) {
            mField = std::vector<std::vector<std::vector<float> > >(
                dimensions,
                std::vector<std::vector<float> >(dimensions, std::vector<float>(dimensions)));
            mIsoFunc = func;
            mIsoVol = isosurfaceVol;
        }

        void Initilise(size_t dimensions,
                       std::function<float(float, float, float)> &&func,
                       float isosurfaceVol) {
            mField = std::vector<std::vector<std::vector<float> > >(
                dimensions,
                std::vector<std::vector<float> >(dimensions, std::vector<float>(dimensions)));
            mIsoFunc = std::move(func);
            mIsoVol = isosurfaceVol;
        }

        void Evaluate(bool useMultithreading = true);

        void Stop() { {
                std::unique_lock<std::mutex> lock(mMutHasFinished);
                mThread.MarkAsFinished();
                mHasFinished.notify_all();
            }
            mThread.Stop();
        }

        bool IsRunning() const {
            return mThread.IsRunning();
        }

        bool HasFinished() const {
            return mThread.HasFinished();
        }

        void Reset() {
            Stop();
            mField.clear();
        }

        std::condition_variable &GetHasFinishedCondition() {
            return mHasFinished;
        }

        std::mutex &GetHasFinishedMutex() {
            return mMutHasFinished;
        }

        const TensorTy &Get() const {
            return mField;
        }

    private:
        float mIsoVol = 0.0f;
        std::function<float(float, float, float)> mIsoFunc;

        std::mutex mMutHasFinished;
        std::condition_variable mHasFinished;

        TensorTy mField;
        DH::Thread mThread;
    };
    //
    //
    //
    class TriSurf {
    public:
        TriSurf() {
        }

        TriSurf(const std::vector<Triangle<float> > &trianglesIn) : mTriangles(trianglesIn) {
        }

        TriSurf(const TriSurf &tri) : mTriangles(tri.mTriangles) {
        }

        TriSurf(TriSurf &&tri) noexcept : mTriangles(std::move(tri.mTriangles)) {
        }

        void Initilise(const std::shared_ptr<ScalarField> &scalarField, float isosurfaceVol, float isosurfaceVal) {
            mField = scalarField;
            mIsoVol = isosurfaceVol;
            mIsoVal = isosurfaceVal;
        }

        void Evaluate(bool useMultithreading = true);

        void Stop() {
            if (mUsingMultithreading)
                mThread.Stop();
        }

        bool IsRunning() const {
            return mThread.IsRunning();
        }

        bool HasFinished() const {
            return mThread.HasFinished();
        }

        void Reset() {
            if (mField != nullptr)
                mField->Stop();

            if (mUsingMultithreading) {
                mThread.Stop();
                mTriangles.clear();
                mField.reset();
            }
        }

        void operator=(const TriSurf &tri) {
            mTriangles = tri.mTriangles;
        }

        void operator=(TriSurf &&tri) noexcept {
            mTriangles = std::move(tri.mTriangles);
        }

        std::vector<Triangle<float> > &GetTriangles() {
            return mTriangles;
        }

        const std::vector<Triangle<float> > &GetTriangles() const {
            return mTriangles;
        }

    private:
        std::shared_ptr<ScalarField> mField;
        std::vector<Triangle<float> > mTriangles;

        float mIsoVol = 0.0f;
        float mIsoVal = 0.0f;

        bool mUsingMultithreading = false;
        DH::Thread mThread;
    };
    //
    //
}
