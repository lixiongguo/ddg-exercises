#pragma once

#include <Eigen/Core>

#include <glm/glm.hpp>

#include <string>
#include <vector>

// 面常向量场的离散 Hodge 分解：
//   ξ = ∇u + n × ∇v + h
// gradient 是可积的梯度部分，curl 是绕法向的旋度部分，harmonic 是两者都吃不下的残余。
struct HodgeSplit {
    std::vector<glm::vec3> gradient;   // 每面 ∇u
    std::vector<glm::vec3> curl;       // 每面 n × ∇v
    std::vector<glm::vec3> harmonic;   // 每面 h
    std::vector<double> potentialU;    // 顶点标量势 u
    std::vector<double> streamV;       // 顶点流函数 v

    double total = 0.0;            // Σ|T| |ξ|²
    double gradientEnergy = 0.0;   // Σ|T| |∇u|²
    double curlEnergy = 0.0;       // Σ|T| |n×∇v|²
    double harmonicEnergy = 0.0;   // Σ|T| |h|²
    double divResidual = 0.0;      // |Div h| / |Div ξ|
    double curlResidual = 0.0;     // |Curl h| / |Curl ξ|

    bool ok = false;
    std::string note;
};

// positions: V×3，faces: F×3（顶点索引），directions: F×3（每面一个切向向量）。
// fold = 4 表示输入是十字场，取分支时会按 90° 的整数倍对齐；普通向量场用 1。
HodgeSplit decomposeHodge(const Eigen::MatrixXd& positions,
                          const Eigen::MatrixXi& faces,
                          const Eigen::MatrixXd& directions,
                          int fold = 4);
