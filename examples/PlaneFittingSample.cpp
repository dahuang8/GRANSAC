#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <random>

#include "GRANSAC.hpp"
#include "PlaneModel.hpp"

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

GRANSAC::VPFloat Slope(int x0, int y0, int x1, int y1)
{
	return (GRANSAC::VPFloat)(y1 - y0) / (x1 - x0);
}

void DrawFullLine(cv::Mat& img, cv::Point a, cv::Point b, cv::Scalar color, int LineWidth)
{
	GRANSAC::VPFloat slope = Slope(a.x, a.y, b.x, b.y);

	cv::Point p(0, 0), q(img.cols, img.rows);

	p.y = -(a.x - p.x) * slope + a.y;
	q.y = -(b.x - q.x) * slope + b.y;

	cv::line(img, p, q, color, LineWidth, cv::LINE_AA, 0);
}

int main(int argc, char * argv[])
{

	double noise = atof(argv[1]);
	double thresh = atof(argv[2]);
	double iter = atoi(argv[3]);
	// create a plane with 100 points
	double a = 1/sqrt(3.0);
	double b = 1/sqrt(3.0);
	double c = 1/sqrt(3.0);
	double d = 50;

	std::cout << "sigma = " << noise << ", threshold = " << thresh << ", iterations = " << iter << std::endl;
	// sample over -10 to 10 over x and y
	std::random_device SeedDevice;
	std::mt19937 RNG = std::mt19937(SeedDevice());
	std::normal_distribution<double> perturb_dist(0, noise);

	std::vector<double> x_row, y_row, z_row;
	std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> CandPoints;
	for (double x = -10.0; x < 10.0; x += 1.0){
		for (double y = -10.0; y < 10.0; y += 1.0) {
			double z = (-d - a * x - b * y) / c;
			// add some noise
			double nx = x + perturb_dist(RNG);
			double ny = y + perturb_dist(RNG);
			double nz = z + perturb_dist(RNG);
			x_row.push_back(nx);
			y_row.push_back(ny);
			z_row.push_back(nz);
			std::shared_ptr<GRANSAC::AbstractParameter> CandPt = std::make_shared<Point3D>(nx,ny,nz);
			CandPoints.push_back(CandPt);
			std::cout << nx << " " << ny << " " << nz << std::endl;
		}
	}

	GRANSAC::RANSAC<Plane3DModel, 3> Estimator;
	Estimator.Initialize(thresh, iter); // Threshold, iterations
	int start = cv::getTickCount();
	Estimator.Estimate(CandPoints);
	int end = cv::getTickCount();
	std::cout << "RANSAC took: " << GRANSAC::VPFloat(end - start) / GRANSAC::VPFloat(cv::getTickFrequency()) * 1000.0 << " ms." << std::endl;

	auto BestInliers = Estimator.GetBestInliers();

	if (BestInliers.size() > 0)
	{
		for (auto& Inlier : BestInliers)
		{
			auto RPt = std::dynamic_pointer_cast<Point3D>(Inlier);
			std::cout << RPt->m_Point3D[0] << " " << RPt->m_Point3D[1] << " " << RPt->m_Point3D[2] << std::endl;
		}
	}

	auto best_plane = Estimator.GetBestModel();
	if (best_plane)
	{
		auto pt0 = std::dynamic_pointer_cast<Point3D>(best_plane->GetModelParams()[0]);
		auto pt1 = std::dynamic_pointer_cast<Point3D>(best_plane->GetModelParams()[1]);
		auto pt2 = std::dynamic_pointer_cast<Point3D>(best_plane->GetModelParams()[2]);
		 // Compute the plane parameters
        double v0[3], v1[3];
        v0[0] = pt1->m_Point3D[0] - pt0->m_Point3D[0];
        v0[1] = pt1->m_Point3D[1] - pt0->m_Point3D[1];
        v0[2] = pt1->m_Point3D[2] - pt0->m_Point3D[2];
        v1[0] = pt2->m_Point3D[0] - pt0->m_Point3D[0];
        v1[1] = pt2->m_Point3D[1] - pt0->m_Point3D[1];
        v1[2] = pt2->m_Point3D[2] - pt0->m_Point3D[2];

        double normal[3];
        normal[0] = v0[1] * v1[2] - v0[2] * v1[1];
        normal[1] = v0[2] * v1[0] - v0[0] * v1[2];
        normal[2] = v0[0] * v1[1] - v0[1] * v1[0];
        double norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        for (int i = 0; i < 3; i++) {
            normal[i] /= norm;
        }

        // ax + by + cz + d = 0
        double m_a = normal[0];
        double m_b = normal[1];
        double m_c = normal[2];
        double m_d = -(m_a * pt0->m_Point3D[0] + m_b * pt0->m_Point3D[1] + m_c * pt0->m_Point3D[2]);
		double ratio = 1;
		std::cout << "orig_a = " << a << ", orig_b = " << b << ", orig_c = " << c << ", orig_d = " << d << std::endl;
		std::cout << "a = " << m_a*ratio << ", b = " << m_b*ratio << ", c = " << m_c*ratio << ", d = " << m_d*ratio << std::endl; 
	}

	/*if (argc != 1 && argc != 3)
	{
		std::cout << "[ USAGE ]: " << argv[0] << " [<Image Size> = 1000] [<nPoints> = 500]" << std::endl;
		return -1;
	}

	int Side = 1000;
	int nPoints = 500;
	if (argc == 3)
	{
		Side = std::atoi(argv[1]);
		nPoints = std::atoi(argv[2]);
	}

	cv::Mat Canvas(Side, Side, CV_8UC3);
	Canvas.setTo(255);

	// Randomly generate points in a 2D plane roughly aligned in a line for testing
	std::random_device SeedDevice;
	std::mt19937 RNG = std::mt19937(SeedDevice());

	std::uniform_int_distribution<int> UniDist(0, Side - 1); // [Incl, Incl]
	int Perturb = 25;
	std::normal_distribution<GRANSAC::VPFloat> PerturbDist(0, Perturb);

	std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> CandPoints;
	for (int i = 0; i < nPoints; ++i)
	{
		int Diag = UniDist(RNG);
		cv::Point Pt(floor(Diag + PerturbDist(RNG)), floor(Diag + PerturbDist(RNG)));
		cv::circle(Canvas, Pt, floor(Side / 100) + 3, cv::Scalar(0, 0, 0), 2, cv::LINE_AA);

		std::shared_ptr<GRANSAC::AbstractParameter> CandPt = std::make_shared<Point2D>(Pt.x, Pt.y);
		CandPoints.push_back(CandPt);
	}

	GRANSAC::RANSAC<Line2DModel, 2> Estimator;
	Estimator.Initialize(20, 100); // Threshold, iterations
	int start = cv::getTickCount();
	Estimator.Estimate(CandPoints);
	int end = cv::getTickCount();
	std::cout << "RANSAC took: " << GRANSAC::VPFloat(end - start) / GRANSAC::VPFloat(cv::getTickFrequency()) * 1000.0 << " ms." << std::endl;

	auto BestInliers = Estimator.GetBestInliers();
	if (BestInliers.size() > 0)
	{
		for (auto& Inlier : BestInliers)
		{
			auto RPt = std::dynamic_pointer_cast<Point2D>(Inlier);
			cv::Point Pt(floor(RPt->m_Point2D[0]), floor(RPt->m_Point2D[1]));
			cv::circle(Canvas, Pt, floor(Side / 100), cv::Scalar(0, 255, 0), -1, cv::LINE_AA);
		}
	}

	auto BestLine = Estimator.GetBestModel();
	if (BestLine)
	{
		auto BestLinePt1 = std::dynamic_pointer_cast<Point2D>(BestLine->GetModelParams()[0]);
		auto BestLinePt2 = std::dynamic_pointer_cast<Point2D>(BestLine->GetModelParams()[1]);
		if (BestLinePt1 && BestLinePt2)
		{
			cv::Point Pt1(BestLinePt1->m_Point2D[0], BestLinePt1->m_Point2D[1]);
			cv::Point Pt2(BestLinePt2->m_Point2D[0], BestLinePt2->m_Point2D[1]);
			DrawFullLine(Canvas, Pt1, Pt2, cv::Scalar(0, 0, 255), 2);
		}
	}

	while (true)
	{
		cv::imshow("RANSAC Example", Canvas);

		char Key = cv::waitKey(1);
		if (Key == 27)
			return 0;
		if (Key == ' ')
			cv::imwrite("LineFitting.png", Canvas);
	}
*/



	return 0;
}
