#pragma once
#include <vector>
#include <cmath>
#include <utility>
//eigen includes
#include <Dense> 

/* NOTES

1. I'm using eigen for linear least squares solving.- Using normal equations because they're intuitive, SVD and QR are better if I need better decimal accuracy but
   for now this is fine. Also I dont understand how SVD or QR work right now so I don't want to implement them if I dont have to
2. For nonlinear stuff, I might have to (probably will have to) do a case by case for each algorithm.
*/

using namespace Eigen;

//regressions
std::vector<float> PolyRegression(std::vector<float> xVals, std::vector<float> yVals, int degree);
std::vector<float> PolyRegressionFixed(std::vector<float> xVals, std::vector<float> yVals, int degree);
//TODO write generalized poly regression for any number of fixed data points
std::vector<std::pair<float, float>> PolyBezierRegression(std::vector<float> xvals, std::vector<float> yvals,int depth, int degree, bool fixed, bool chord_length);

//helper functions for bezier regressions
std::vector<std::pair<float, float>>  ComputeBez(std::vector<float> xvals, std::vector<float> yvals, VectorXf tvals, int degree, bool fixed);
VectorXf Reparameterize(std::vector<float> xvals, std::vector<float> yvals, std::vector<std::pair<float, float>> cpoints, VectorXf tvals);
float NewtonRootFinder(float xval, float yval, std::vector<std::pair<float, float>> cpoints, float t_index, float tupper, float tlower);
std::pair<float, float> BezValue(std::vector<std::pair<float, float>> cpoints, float t_index, int degree);

//helper math functions
int fact(int num);
int Choose(int k, int n);


//non matrix implementation of cubic bezier 
//OLD, BAD, DONT USE
template <typename T, typename K>
std::vector<std::pair<float, float>> CubicBezierTwo(std::vector<T> xVals, std::vector<K> yVals, bool fixed = true) {
	std::vector<std::pair<float, float>> control_points;
	float A1 = 0;
	float A2 = 0;
	float A12 = 0;
	float C1X = 0;
	float C1Y = 0;
	float C2X = 0;
	float C2Y = 0;
	float t = 0;
	for (int i = 0; i < xVals.size(); i++) {
		t = static_cast<float>(i * 1.0 / xVals.size());
		A1 = A1 + 9 * pow(t, 2) * pow(1 - t, 4);
		A2 = A2 + 9 * pow(t, 4) * pow(1 - t, 2);
		A12 = A12 + 9 * pow(t, 3) * pow(1 - t, 3);
		C1X = C1X + 3 * t * pow(1 - t, 2) * (xVals[i] - xVals[0] * pow(1 - t, 3) - xVals.back() * pow(t, 3));
		C1Y = C1Y + 3 * t * pow(1 - t, 2) * (yVals[i] - yVals[0] * pow(1 - t, 3) - yVals.back() * pow(t, 3));
		C2X = C2X + 3 * pow(t, 2) * (1 - t) * (xVals[i] - xVals[0] * pow(1 - t, 3) - xVals.back() * pow(t, 3));
		C2Y = C2Y + 3 * pow(t, 2) * (1 - t) * (yVals[i] - yVals[0] * pow(1 - t, 3) - yVals.back() * pow(t, 3));
	}
	float denom = A1 * A2 - pow(A12, 2);
	if (denom == 0) {
		denom = 1;
	}


	control_points.push_back(std::make_pair((A2 * C1X - A12 * C2X) / denom, (A2 * C1Y - A12 * C2Y) / denom));
	control_points.push_back(std::make_pair((A1 * C2X - A12 * C1X) / denom, (A1 * C2Y - A12 * C1Y) / denom));
	return control_points;
}

