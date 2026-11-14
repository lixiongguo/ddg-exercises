#pragma once

#include <string>

namespace polyscope {
class SurfaceMesh;
class SurfaceVectorQuantity;
}

enum FieldStep {
    StepCurvature = 0,
    StepPrincipal = 1,
    StepRosy = 2,
    StepPoly = 3,
    StepTree = 4,
    StepGenerators = 5,
    StepOmega = 6,
    StepExact = 7,
    StepCoexact = 8,
    StepHarmonic = 9,
    StepBasis = 10,
    StepConnection = 11,
    StepCount = 12
};

class FieldProcess {
  public:
    bool prepare(const std::string& path, std::string& error);
    void newOneForm();
    void nextBasis();
    void resetSingularities();
    bool computeConnection(std::string& error);
    void suspendPicking();
    void showSelected();
    void draw(int process, bool resetCamera, std::string& status);

    int euler() const { return euler_; }
    double singularitySum() const;
    double* singularityAt(int vertex);
    int basisIndex() const { return basisIndex_; }
    int basisCount() const;
    polyscope::SurfaceVectorQuantity* connectionForm() const { return connectionForm_; }
    bool& showConnectionForm() { return showConnectionForm_; }

  private:
    void decompose();
    void drawDirection(polyscope::SurfaceMesh* mesh);

    bool ready_ = false;
    bool picking_ = false;
    bool showConnectionForm_ = false;
    int euler_ = 0;
    int basisIndex_ = 0;
    double lengthScale_ = 1.0;
    std::string path_;
    std::string note_;
    polyscope::SurfaceMesh* drawn_ = nullptr;
    polyscope::SurfaceVectorQuantity* connectionForm_ = nullptr;

    struct State;
    State* state_ = nullptr;
};

FieldProcess& fieldProcess();
