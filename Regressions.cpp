#include "Regressions.h"

/* REGRESSION LIBRARY
*
* This library is designed to help people do linear regressions on sets of data
* I wrote this as part of another project, and it turned into its own thing.
* I try to add as much flexibility as possible without making it super
* convoluted, hopefully people can get value from this. Otherwise, at least I did
*
* RETURN VALUES: Everything is returned as a vector of floats, or vector of pairs for parametric equations
*
* DEPENDENCIES: I use eigen to optimize the math using matrix multiplication.
* Not sure which versions are compatible, sorry. I may try to write alternatives that
*  are dependancy free, but then they will be slow cuz IDK how to optimize matrix math
*
*
*/


/* DESCRIPTION: Linear regression on a set of data. Input your data and the degree polynomial you want,
* the function should take care of the rest :)
*
* NOTE: Due to (I think) floating point imprecision, currently this only works up to polynomial degree 6
*
* PARAMS:
* vector<float> xvals and yvals are the input data, xvals are the regressor and yvals are the regressand
* int degree is the degree polynomial, with degree 0 being y = a, degree 1 being y = ax + b, and so on
*
* RETURN VALUE: this function returns a vector of the coefficients, with the length == degree + 1
*/
std::vector<float> PolyRegression(std::vector<float> xVals, std::vector<float> yVals, int degree) {
	std::vector<float> coefficientsVec;
	int n = xVals.size();
	int m = degree;

	assert(n > m, "data has to be overdetermined");
	assert(m >= 0, "regression degree must be positive");
	assert(n == yVals.size(), "data columns must be of same length");

	//generating x Matrix
	MatrixXf   xMat(xVals.size(), degree + 1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= m; j++) {
			xMat(i, j) = pow(xVals[i], j);
		}
	}
	//generating y Matrix
	VectorXf yVec(n);
	for (int j = 0; j < n; j++)
		yVec(j) = yVals[j];


	//solving for least squares. Uses equation derived here https://en.wikipedia.org/wiki/Proofs_involving_ordinary_least_squares
	VectorXf coefficients = (xMat.transpose() * xMat).ldlt().solve(xMat.transpose() * yVec);
	for (float val : coefficients)
		coefficientsVec.push_back(val);
	return coefficientsVec;
}


//Poly regression with fixed endpoints
//TODO make this generalized to work with any number of point constraints
//TODO: fix documetnation for this function
std::vector<float> PolyRegressionFixed(std::vector<float> xVals, std::vector<float> yVals, int degree) {
	std::vector<float> coefficientsVec;
	int n = xVals.size();
	int m = degree;

	assert(n > m, "data has to be overdetermined");
	assert(m >= 0, "regression degree must be positive");
	assert(n == yVals.size(), "data columns must be of same length");
	assert(n >= 2, "must have at least two points for fixed endpoint linreg");
	float x0 = xVals[0];
	float xf = xVals.back();
	float y0 = yVals[0];
	float yf = yVals.back();

	//generating x Matrix
	MatrixXf   xMat(xVals.size(), degree + 1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= m; j++) {
			xMat(i, j) = pow(xVals[i], j);
		}
	}
	//generating y Matrix
	VectorXf yVec(n);
	for (int j = 0; j < n; j++) {
		float inval = yVals[j];
		yVec(j) = inval;
	}



	//resultant matrix
	VectorXf resultant(m + 1);
	resultant = xMat.transpose() * yVec;
	resultant.conservativeResize(resultant.rows() + 2, resultant.cols());
	resultant(m + 1) = y0;
	resultant(m + 2) = yf;

	//solving for least squares. Uses equation  derived here http://www.seas.ucla.edu/~vandenbe/133A/lectures/cls.pdf 
	xMat = xMat.transpose() * xMat;

	xMat.conservativeResize(xMat.rows() + 2, xMat.cols() + 2);

	//generating constraint vectors
	VectorXf beg_column = VectorXf::Zero(xMat.rows());
	VectorXf end_column = VectorXf::Zero(xMat.rows());

	for (int i = 0; i < m + 1; i++) {
		beg_column(i) = pow(x0, i);
		end_column(i) = pow(xf, i);
	}
	std::vector<float> beg_test;
	std::vector<float> end_test;
	std::vector<float> result_test;
	std::vector<std::vector<float>> mat_test;
	std::vector<float> temp_mat;

	//im pretty sure for the columns im supposed to divide by 2
	//but i tested it and i got zero difference in coefficient values
	//Either floating point imprecision means it doesnt matter, or my math is wrong
	//Regardless, the results look good... may look into it more in the future when its not 3 AM
	xMat.col(xMat.cols() - 2) = beg_column;
	xMat.col(xMat.cols() - 1) = end_column;
	xMat.row(xMat.rows() - 2) = beg_column.transpose();
	xMat.row(xMat.rows() - 1) = end_column.transpose();
	for (int i = 0; i < xMat.rows(); i++) {
		beg_test.push_back(beg_column(i));
		end_test.push_back(end_column(i));
		mat_test.push_back(temp_mat);
		for (int j = 0; j < xMat.cols(); j++) {
			mat_test[i].push_back(xMat(i, j));
		}

	}
	VectorXf coefficients = xMat.ldlt().solve(resultant);
	for (float val : coefficients)
		coefficientsVec.push_back(val);
	//get rid of the two z values
	coefficientsVec.pop_back();
	coefficientsVec.pop_back();
	return coefficientsVec;
}

/* DESCRIPTION: poly bezier regression Fits a set of data points to a bezier curve using a linear regression matrix
 * The regression goes back and forth between fitting and conditioning the regressor variable t
 *
 NOTE: Due to (I think) floating point imprecision, currently this only works up to degree 11 consistently
 *
 * PARAMS:
 * vector<float> xvals and yvals are the input vectors describing the data to be regressed upon
 * int depth is the number of reparameterizations done to condition the t values for newton approximation
 * int degree is the degree of the bezier curve (degree 3 -> cubic)
 * bool fixed is whether or not the endpoints need to stay fixed
 * bool chord_length is whether or not the t values should be originally chord-length spaced or not (alternative is linear spacing)
 *
 * RETURN VALUE: The function returns a vector of pairs which are the bezier control points
*/
std::vector<std::pair<float, float>> PolyBezierRegression(std::vector<float> xvals, std::vector<float> yvals, int depth, int degree, bool fixed, bool chord_length) {
	int xsize = xvals.size();
	int ysize = yvals.size();
	std::vector<std::pair<float, float>> control_points;

	//assertions so program doesnt break
	assert(xsize == ysize, "vectors must be the same size");
	assert(xsize >= 2, "must have at least 2 data points");
	assert(degree > 0, "degree must be greater than 0");
	assert(depth >= 0, "depth must be non-negative");

	//for 2 points, put middle points on end points to form a straight line
	if (xsize == 2) {
		control_points.push_back({ xvals[0],yvals[0] });
		control_points.push_back({ xvals[0],yvals[0] });
		control_points.push_back({ xvals[1],yvals[1] });
		control_points.push_back({ xvals[1],yvals[1] });
		return control_points;
	}

	//creating original regressor matrix
	VectorXf tvals;
	std::vector<float> path_length;
	if (!chord_length) {
		tvals = VectorXf::LinSpaced(xsize, 0, 1);
	}
	else {
		tvals = VectorXf(xsize);
		path_length.push_back(0);
		for (int i = 1; i < xsize; i++) {
			float distance = sqrt(pow(xvals[i] - xvals[i - 1], 2) + pow(yvals[i] - yvals[i - 1], 2));
			path_length.push_back(path_length.back() + distance);
		}
		for (int i = 0; i < xsize; i++) {
			float tLen = path_length[i] / path_length.back();
			tvals(i) = tLen;
		}
	}

	//iterations for reparameterization
	if (degree >= 2)
		for (int i = 0; i < depth; i++) {
			control_points = ComputeBez(xvals, yvals, tvals, degree, fixed); //compute bezier and store control points
			tvals = Reparameterize(xvals, yvals, control_points, tvals); //reparamaterize and store T values
		}
	control_points = ComputeBez(xvals, yvals, tvals, degree, fixed);
	return control_points;
}

/* DESCRIPTION: ComputeBez is a helper function that does the control point computation for the PolyBezier function
 *
 * PARAMS:
 * vector<float> xvals and yvals are the input vectors describing the data to be regressed upon
 * VectorXf tvals is an eigen datatype, and it stores the t values for the bezier
 * bool degree is the degree of the bezier
 * bool fixed is whether or not the endpoints need to stay fixed
 *
 * RETURN VALUE: The function returns a vector of pairs which are the control points
*/
std::vector<std::pair<float, float>> ComputeBez(std::vector<float> xvals, std::vector<float> yvals, VectorXf tvals, int degree, bool fixed) {
	float p0_x = xvals.front();
	float p3_x = xvals.back();
	float p0_y = yvals.front();
	float p3_y = yvals.back();
	int xsize = xvals.size();
	MatrixXf var_matrix;
	VectorXf x_controls;
	VectorXf y_controls;
	std::vector<std::pair<float, float>> control_points;

	if (fixed) {
		//creating regressor matrix
		MatrixXf regressors(xsize, degree - 1);
		for (int i = 0; i < degree - 1; i++) {
			float choosey = Choose(degree, i + 1);
			for (int j = 0; j < xsize; j++)
			{
				regressors(j, i) = pow(1 - tvals(j), degree - 1 - i) * pow(tvals(j), i + 1) * choosey;
			}
		}
		//creating offset regressands

		VectorXf y_vec(xsize);
		for (int j = 0; j < xsize; j++)
			y_vec(j) = yvals[j] - p0_y * pow(1 - tvals(j), degree) - p3_y * pow(tvals(j), degree);

		VectorXf x_vec(xsize);
		for (int j = 0; j < xsize; j++)
			x_vec(j) = xvals[j] - p0_x * pow(1 - tvals(j), degree) - p3_x * pow(tvals(j), degree);

		//performing calculations to find the control points
		var_matrix = ((regressors.transpose() * regressors).inverse()) * regressors.transpose();
		x_controls = var_matrix * x_vec;
		y_controls = var_matrix * y_vec;

		control_points.push_back(std::make_pair(xvals.front(), yvals.front()));
		for (int i = 0; i < degree - 1; i++) {
			control_points.push_back(std::make_pair(x_controls(i), y_controls(i)));
		}
		control_points.push_back(std::make_pair(xvals.back(), yvals.back()));
	}
	else {
		//same format as the fixed endpoints but including the endpoints as coefficients
		MatrixXf regressors(xsize, degree + 1);
		for (int i = 0; i < degree + 1; i++) {
			float choosey = Choose(degree, i);
			for (int j = 0; j < xsize; j++)
			{
				regressors(j, i) = pow(1 - tvals(j), degree - i) * pow(tvals(j), i) * choosey;
			}
		}
		//no offset because the endpoints are also coefficients
		VectorXf y_vec(xsize);
		for (int j = 0; j < xsize; j++)
			y_vec(j) = yvals[j];

		VectorXf x_vec(xsize);
		for (int j = 0; j < xsize; j++)
			x_vec(j) = xvals[j];

		var_matrix = ((regressors.transpose() * regressors).inverse()) * regressors.transpose();
		x_controls = var_matrix * x_vec;
		y_controls = var_matrix * y_vec;
		for (int i = 0; i < degree + 1; i++) {
			control_points.push_back(std::make_pair(x_controls(i), y_controls(i)));
		}
	}
	return control_points;

}


/* DESCRIPTION: Reparameterize is a helper function that does the control point computation for the PolyBezier function
 *
 * PARAMS:
 * vector<float> xvals and yvals are the input vectors describing the data to be regressed upon
 * vector<pair<float,float>> cpoints is the vector of control points for a given bezier curve
 * VectorXf tvals is an eigen datatype, and it stores the t values for the bezier
 *
 * RETURN VALUE: The function returns a reparameterized VectorXf of tvals
*/
VectorXf Reparameterize(std::vector<float> xvals, std::vector<float> yvals, std::vector<std::pair<float, float>> cpoints, VectorXf tvals) {
	int size = xvals.size();
	VectorXf reparam_t(size);
	std::vector<float> t_test;
	reparam_t(0) = 0.000001; //to prevent floating point error i have to make this positive
	t_test.push_back(reparam_t(0));

	for (int i = 1; i < size - 1; i++) {
		reparam_t(i) = NewtonRootFinder(xvals[i], yvals[i], cpoints, tvals(i), tvals(i + 1), tvals(i - 1));
		t_test.push_back(reparam_t(i));
	}
	reparam_t(size - 1) = tvals(size - 1);
	return reparam_t;
}

/* DESCRIPTION: NewtonRootFinder is a function that approximates roots. This specific function
 * is set up only to approximate roots for bezier curves. Trying to generalize this to any curve
 * would be challenging, so currently it can only approximate roots for bezier curves
 *
 * PARAMS:
 * vector<float> xvals and yvals are the input vectors describing the data to be regressed upon
 * vector<pair<float,float>> cpoints is the vector of control points for a given bezier curve
 * float t_index is the current t value to be reparameterized
 * float tupper is the next t value in the curve, used as a ceiling for t_index
 * float tlower is the previous t value in the curve, used as a floor for t_index
 *
 * RETURN VALUE: The function returns a reparameterized tval
*/
float NewtonRootFinder(float xval, float yval, std::vector<std::pair<float, float>> cpoints, float t_index, float tupper, float tlower) {
	int degree = cpoints.size() - 1;
	std::pair<float, float> base = BezValue(cpoints, t_index, degree);
	std::vector<std::pair<float, float>> prime_control;
	//calculating derivatives
	for (int i = 0; i <= degree - 1; i++) {
		prime_control.push_back({ degree * (cpoints[i + 1].first - cpoints[i].first), degree * (cpoints[i + 1].second - cpoints[i].second) });
	}
	std::vector<std::pair<float, float>> dbl_prime_control;
	for (int i = 0; i <= degree - 2; i++) {
		dbl_prime_control.push_back({ (degree - 1) * (prime_control[i + 1].first - prime_control[i].first),  (degree - 1) * (prime_control[i + 1].second - prime_control[i].second) });
	}
	std::pair<float, float> prime = BezValue(prime_control, t_index, degree - 1);
	std::pair<float, float> dbl_prime = BezValue(dbl_prime_control, t_index, degree - 2);

	// new t = t - f(x)/f'(x)
	//our f(x) is the derivative of the squared distance from the point to the bezier curve
	float numerator = (base.first - xval) * prime.first + (base.second - yval) * prime.second;
	float denominator = pow(prime.first, 2) + pow(prime.second, 2) + (base.first - xval) * dbl_prime.first + (base.second - yval) * dbl_prime.second;
	float fraction = numerator / denominator;
	if (denominator == 0)
		return t_index;
	//important to clamp t values such that they dont cross each other
	if (t_index - fraction < tlower)
		t_index = tlower;
	else if (t_index - fraction > tupper)
		t_index = tupper;
	else
		t_index = t_index - fraction;
	return t_index;
}

/* DESCRIPTION: BezValue is a function that calculates the x and y values of a point on a bezier curve
 * at a given time t
 *
 * PARAMS:
 * vector<pair<float,float>> cpoints is the vector of control points for a given bezier curve
 * float t_index is the current t value used to calculate the x and y values
 * int degree is the degree of the bezier curve
 *
 * RETURN VALUE: The function returns the x and y coordinates of a bezier curve with cpoints at time t_index
*/
std::pair<float, float> BezValue(std::vector<std::pair<float, float>> cpoints, float t_index, int degree) {
	std::pair<float, float> bez_value = std::make_pair(0, 0);
	for (int i = 0; i <= degree; i++) {
		bez_value.first += Choose(degree, i) * pow(t_index, i) * pow(1 - t_index, degree - i) * cpoints[i].first;
		bez_value.second += Choose(degree, i) * pow(t_index, i) * pow(1 - t_index, degree - i) * cpoints[i].second;
	}
	return bez_value;
}


/* DESCRIPTION: fact is a helper math function that calculates the factorial of an integer num
 *
 * PARAMS:
 * int num is the num to be factorialized
 *
 * RETURN VALUE: integer that is the factorial of num
*/
int fact(int num) {
	for (int i = num - 1; i > 0; i--) {
		num = num * i;
	}
	if (num == 0)
		num = 1;
	return num;
}


/* DESCRIPTION: Choose is an implementation of the kCn operation, i.e. combinations
 *
 * PARAMS:
 * int k is the number that chooses
 * int n is the number that is chosen
 *
 * RETURN VALUE: int that is the result of performing kCn
*/
int Choose(int k, int n) {
	return (fact(k) / (fact(n) * fact(k - n)));
}






