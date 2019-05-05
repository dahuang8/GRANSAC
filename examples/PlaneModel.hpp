#pragma once

#include "AbstractModel.hpp"

typedef std::array<GRANSAC::VPFloat, 3> Vector3VP;

class Point3D : public GRANSAC::AbstractParameter {
 public:
    Point3D(GRANSAC::VPFloat x, GRANSAC::VPFloat y, GRANSAC::VPFloat z) {
        m_Point3D[0] = x;
        m_Point3D[1] = y;
        m_Point3D[2] = z;
    };

    Vector3VP m_Point3D;
};

class Plane3DModel : public GRANSAC::AbstractModel<3> {
 protected:
    // Parametric form
    GRANSAC::VPFloat m_a, m_b, m_c, m_d;  // ax + by + cz + d = 0
    GRANSAC::VPFloat m_DistDenominator;   // = sqrt(a^2 + b^2 + c^2). Stored for efficiency reasons

    virtual GRANSAC::VPFloat ComputeDistanceMeasure(std::shared_ptr<GRANSAC::AbstractParameter> Param) override {
        auto ExtPoint3D = std::dynamic_pointer_cast<Point3D>(Param);
        if (ExtPoint3D == nullptr)
            throw std::runtime_error(
                "Plane3DModel::ComputeDistanceMeasure() - Passed "
                "parameter are not of type Point3D.");

        // Return distance between passed "point" and this line
        // http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
        GRANSAC::VPFloat Numer = fabs(m_a * ExtPoint3D->m_Point3D[0] + m_b * ExtPoint3D->m_Point3D[1] +
                                      m_c * ExtPoint3D->m_Point3D[2] + m_d);
        GRANSAC::VPFloat Dist = Numer / m_DistDenominator;

        //// Debug
        // std::cout << "Point: " << ExtPoint3D->m_Point2D[0] << ", " <<
        // ExtPoint3D->m_Point2D[1] << std::endl;
        // std::cout << "Line: " << m_a << " x + " << m_b << " y + "  << m_c <<
        // std::endl;
        // std::cout << "Distance: " << Dist << std::endl << std::endl;

        return Dist;
    };

 public:
    Plane3DModel(const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> &InputParams) {
        Initialize(InputParams);
    };

    virtual void Initialize(const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> &InputParams) override {
        if (InputParams.size() != 3)
            throw std::runtime_error(
                "Plane3DModel - Number of input parameters does "
                "not match minimum number required for this "
                "model.");

        // Check for AbstractParamter types
        auto Point1 = std::dynamic_pointer_cast<Point3D>(InputParams[0]);
        auto Point2 = std::dynamic_pointer_cast<Point3D>(InputParams[1]);
        auto Point3 = std::dynamic_pointer_cast<Point3D>(InputParams[2]);
        if (Point1 == nullptr || Point2 == nullptr || Point3 == nullptr)
            throw std::runtime_error("Line2DModel - InputParams type mismatch. It is not a Point3D.");

        std::copy(InputParams.begin(), InputParams.end(), m_MinModelParams.begin());

        // Compute the plane parameters
        double v0[3], v1[3];
        v0[0] = Point2->m_Point3D[0] - Point1->m_Point3D[0];
        v0[1] = Point2->m_Point3D[1] - Point1->m_Point3D[1];
        v0[2] = Point2->m_Point3D[2] - Point1->m_Point3D[2];
        v1[0] = Point3->m_Point3D[0] - Point1->m_Point3D[0];
        v1[1] = Point3->m_Point3D[1] - Point1->m_Point3D[1];
        v1[2] = Point3->m_Point3D[2] - Point1->m_Point3D[2];

        double normal[3];
        normal[0] = v0[1] * v1[2] - v0[2] * v1[1];
        normal[1] = v0[2] * v1[0] - v0[0] * v1[2];
        normal[2] = v0[0] * v1[1] - v0[1] * v1[0];
        double norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        for (int i = 0; i < 3; i++) {
            normal[i] /= norm;
        }

        // mx - y + d = 0
        m_a = normal[0];
        m_b = normal[1];
        m_c = normal[2];
        m_d = -(m_a * Point1->m_Point3D[0] + m_b * Point1->m_Point3D[1] + m_c * Point1->m_Point3D[2]);

        m_DistDenominator = sqrt(m_a * m_a + m_b * m_b + m_c * m_c);  // Cache square root for efficiency
    };

    virtual std::pair<GRANSAC::VPFloat, std::vector<std::shared_ptr<GRANSAC::AbstractParameter>>> Evaluate(
        const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> &EvaluateParams, GRANSAC::VPFloat Threshold) {
        std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> Inliers;
        int nTotalParams = EvaluateParams.size();
        int nInliers = 0;

        for (auto &Param : EvaluateParams) {
            if (ComputeDistanceMeasure(Param) < Threshold) {
                Inliers.push_back(Param);
                nInliers++;
            }
        }

        GRANSAC::VPFloat InlierFraction =
            GRANSAC::VPFloat(nInliers) / GRANSAC::VPFloat(nTotalParams);  // This is the inlier fraction

        return std::make_pair(InlierFraction, Inliers);
    };
};
