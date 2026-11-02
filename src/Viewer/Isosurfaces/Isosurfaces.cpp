//Adapted from https://github.com/nihaljn/marching-cubes

#include "Isosurfaces.h"
#include <chrono>
#include "Dahlia/Tools/Threading.h"

//
//
namespace SCESS::Viewer {
    //
    //
    using namespace DH;
    //
    //
#pragma region Marching Cubes
    //
    int MarchingCubes::CalculateCubeIndex(GridCell &cell, float isovalue) {
        int cubeIndex = 0;
        for (int i = 0; i < 8; i++)
            if (cell.value[i] < isovalue) cubeIndex |= (1 << i);
        return cubeIndex;
    }
    //
    void MarchingCubes::GenerateIntersectionCoordinates(GridCell &cell, float isovalue) {
        int cubeIndex = CalculateCubeIndex(cell, isovalue);
        int intersectionsKey = edgeTable[cubeIndex];

        int idx = 0;
        while (intersectionsKey) {
            if (intersectionsKey & 1) {
                int v1 = edgeToVertices[idx].first, v2 = edgeToVertices[idx].second;
                glm::vec3 intersectionPoint = Interpolate(cell.vertex[v1], cell.value[v1],
                                                          cell.vertex[v2], cell.value[v2], isovalue);
                mIntersections[idx] = intersectionPoint;
            }
            idx++;
            intersectionsKey >>= 1;
        }
    }
    //
    glm::vec3 MarchingCubes::Interpolate(glm::vec3 &v1, float val1, glm::vec3 &v2, float val2, float isovalue) {
        glm::vec3 interpolated;
        float mu = (isovalue - val1) / (val2 - val1);

        interpolated.x = mu * (v2.x - v1.x) + v1.x;
        interpolated.y = mu * (v2.y - v1.y) + v1.y;
        interpolated.z = mu * (v2.z - v1.z) + v1.z;

        return interpolated;
    }
    //
    void MarchingCubes::TriangulateIntersections(std::vector<Triangle<float> > &triangles, int cubeIndex) {
        for (int i = 0; triangleTable[cubeIndex][i] != -1; i += 3) {
            Triangle<float> triangle;
            for (int j = 0; j < 3; j++)
                triangle[j] = mIntersections[triangleTable[cubeIndex][i + j]];
            triangles.push_back(triangle);
        }
    }
    //
    void MarchingCubes::TriangulateCell(std::vector<Triangle<float> > &trianglesOut,
                                        GridCell &cell, float isovalue) {
        int cubeIndex = CalculateCubeIndex(cell, isovalue);
        GenerateIntersectionCoordinates(cell, isovalue);
        TriangulateIntersections(trianglesOut, cubeIndex);
    }
    //
    void MarchingCubes::TriangulateField(
        const ScalarField &field,
        float isovalue, std::vector<Triangle<float> > &trianglesOut, const Thread *thread) {
        const ScalarField::TensorTy &scalarF = field.Get();
        size_t max = scalarF.size();

        for (size_t i = 0; i + 1 < max; i++) {
            for (size_t j = 0; j + 1 < max; j++) {
                if (thread != nullptr)
                    if (thread->ShouldStop())
                        return;
                for (size_t k = 0; k + 1 < max; k++) {
                    float x = i, y = j, z = k;
                    // cell ordered according to convention in referenced website
                    GridCell cell =
                    {
                        {
                            {x, y, z}, {x + 1.0f, y, z},
                            {x + 1.0f, y, z + 1.0f}, {x, y, z + 1.0f},
                            {x, y + 1.0f, z}, {x + 1.0f, y + 1.0f, z},
                            {x + 1.0f, y + 1.0f, z + 1.0f}, {x, y + 1.0f, z + 1.0f}
                        },
                        {
                            scalarF[i][j][k], scalarF[i + 1][j][k],
                            scalarF[i + 1][j][k + 1], scalarF[i][j][k + 1],
                            scalarF[i][j + 1][k], scalarF[i + 1][j + 1][k],
                            scalarF[i + 1][j + 1][k + 1], scalarF[i][j + 1][k + 1]
                        }
                    };

                    TriangulateCell(trianglesOut, cell, isovalue);
                }
            }
        }
    }
    //
#pragma endregion
    //
#pragma region Isosurface
    //
    void ScalarField::Evaluate(bool useMultithreading) {
        {
            std::lock_guard<std::mutex> lock(mMutHasFinished);
            mEvaluationDone = false;
            mEvaluationCompleted = false;
        }

        auto workerFunc = [this](bool useMultithreading) {
            size_t resolution = mField.size();
            bool completed = true;

            for (size_t i = 0; i < resolution && completed; i++) {
                if (useMultithreading && mThread.ShouldStop()) {
                    completed = false;
                    break;
                }
                for (size_t j = 0; j < resolution; j++) {
                    for (size_t k = 0; k < resolution; k++) {
                        glm::vec3 pos(
                            (2.0 * float(i) / float(resolution) - 1.0) * mIsoVol,
                            (2.0 * float(j) / float(resolution) - 1.0) * mIsoVol,
                            (2.0 * float(k) / float(resolution) - 1.0) * mIsoVol);
                        mField[i][j][k] = mIsoFunc(pos.x, pos.y, pos.z);
                    }
                }
            }

            {
                std::lock_guard<std::mutex> lock(mMutHasFinished);
                mEvaluationDone = true;
                mEvaluationCompleted = completed;
            }
            mHasFinished.notify_all();

            if (useMultithreading)
                mThread.MarkAsFinished();
        };

        if (!useMultithreading) {
            workerFunc(false);
            return;
        }

        mThread.Start(workerFunc, useMultithreading);
    }
    //
    bool ScalarField::WaitUntilEvaluated(const DH::Thread *waiter) {
        using namespace std::chrono_literals;
        std::unique_lock<std::mutex> lock(mMutHasFinished);

        //Wait on the evaluation state rather than the notification alone so a notification sent before waiting is
        //not missed, and poll the waiter so that it can be stopped while waiting
        while (!mEvaluationDone) {
            if (waiter && waiter->ShouldStop())
                return false;
            mHasFinished.wait_for(lock, 20ms);
        }

        return mEvaluationCompleted;
    }
    //
    void TriSurf::Evaluate(bool useMultithreading) {
        auto workerFunc = [this]() {
            std::unique_ptr<MarchingCubes> marchingCubes = std::make_unique<MarchingCubes>();

            mTriangles.clear();

            if (mUsingMultithreading && !mField->WaitUntilEvaluated(&mThread)) {
                mThread.MarkAsFinished();
                return;
            }

            marchingCubes->TriangulateField(*mField, mIsoVal, mTriangles, mUsingMultithreading ? &mThread : nullptr);

            //Rescale vertices
            for (size_t i = 0; i < mTriangles.size(); i++) {
                for (uint8_t j = 0; j < 3; j++) {
                    mTriangles[i][j] = ((mTriangles[i][j] * 2.0f / float(mField->Get().size())) - 1.0f) * mIsoVol;
                }
            }

            if (mUsingMultithreading) mThread.MarkAsFinished();
        };

        mUsingMultithreading = useMultithreading;
        if (!mUsingMultithreading) {
            workerFunc();
            return;
        }

        mThread.Start(workerFunc);
    }
    //
#pragma endregion
    //
    //
}
