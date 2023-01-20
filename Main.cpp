#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <algorithm>
#include <iostream>
#include <functional>
#include <map>
#include <GLFW/glfw3.h>
#include <vector>
#include <stack>
#include <queue>
#include "Regressions.h"

using namespace std;

//forward declarations
float GetSlope(ImVec2 p1, ImVec2 p2);

//enums that describe fragment data
//These enums are not set up to be typesafe, just a pain to implement that and also reduces code readability
//Since these are all internal anyways it should be fine unless I make a mistake (which is likely tbf)
enum Regression
{
	None = 1 << 0,
	FlatLine = 1 << 1,
	VertLine = 1 << 2,
	Linear = 1 << 3,
	Quadratic = 1 << 4,
	PolyFit = 1 << 5,
	BezLinear = 1 << 6,
	BezQuadratic = 1 << 7,
	BezCubic = 1 << 8,
	BezPolyFit = 1 << 9
};
enum Direction
{
	Forward = 1 << 0,
	Backward = 1 << 1,
	Mixed = 1 << 2
};
enum WD_Window
{
	Drawing = 1 << 0,
	Shaping = 1 << 1,
	Superposition = 1 << 2,
	Playback = 1 << 3
};

enum DivType
{
	Manual = 1 << 0,
	Extrema = 1 << 1,
	ChangeDir = 1 << 2,
	Other = 1 << 3
};



//enum operator overloads for comparison
inline DivType operator&(DivType a, DivType b) {
	return static_cast<DivType>(static_cast<int>(a) & static_cast<int>(b));
}
inline DivType operator|(DivType a, DivType b) {
	return static_cast<DivType>(static_cast<int>(a) | static_cast<int>(b));
}

//ImVec2 equality operator overload (arbitrary min difference for floating point imprecision correction)
bool operator ==(const ImVec2& lhs, const ImVec2& rhs) {
	float diff = fabs(lhs.x - rhs.x) + fabs(lhs.y - rhs.y);
	float lhs_max = fabs(lhs.x + lhs.y);
	float rhs_max = fabs(rhs.x + rhs.y);
	float max_val = max(rhs_max, lhs_max);

	if (diff <= max_val * FLT_EPSILON * 10)
		return true;
	return false;
}

//ImVec2 addition overload
ImVec2 operator +(const ImVec2& lhs, const ImVec2& rhs) {
	return ImVec2(lhs.x + rhs.x, lhs.y + rhs.y);
}

//ImVec2 subtraction overload
ImVec2 operator -(const ImVec2& lhs, const ImVec2& rhs) {
	return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y);
}

//  struct for individual line segments to save state and allow for manipulation
//  TODO: add curated drawing options, like straight line drawing, bezier drawing,
//  spline drawing, etc. These drawing modes will have the raw mouse input values
//  be regressed, or straightened into a line, real time, and then the points will
//  be saved into a draw fragment at the end of each frame post manipulation
//  that way, raw_data represents the raw visual data, not technically raw mouse data
//  this also prevents the need for a million constructors

//  a note on RAW_DATA: raw_data is essentially a reset state for the *points* data. It is not used in any comparisons (except for == operator overload)
//  or any functions that assess data, because it should never be used in the program, except as a reset state/backup. Thus it is the responsibility of
//  the rest of the code (the main method mainly) to be aware of this fact, and NEVER do anything with raw data other than save/load it
class DrawFragment {
public:
	vector<ImVec2> points;   // points used for whats displayed on screen
	vector<ImVec2> raw_data; // raw data of what's drawn on the page, used as a backup
	int reg_degree;          // regression degree
	Regression reg_type;
	Direction dir;
	DivType div;
	DrawFragment(vector<ImVec2> points, Direction dir, DivType div, Regression reg_type, int reg_degree) {
		this->div = div;
		this->points = points;
		this->raw_data = points;
		this->reg_type = reg_type;
		this->reg_degree = reg_degree;
		this->dir = dir;
	}
	DrawFragment(vector<ImVec2> points, Direction dir, DivType div) {
		this->points = points;
		this->raw_data = points;
		this->dir = dir;
		this->reg_type = None;
		this->reg_degree = -1;
		this->div = div;
	}
	DrawFragment(vector<ImVec2> points) {
		this->points = points;
		this->raw_data = points;
		this->reg_type = None;
		this->reg_degree = -1;
		this->dir = Mixed;
		this->div = Other;
	}
	DrawFragment() {}
	~DrawFragment() {}

	vector<float> Xvals() {
		vector<float> vec;
		for (ImVec2 p : this->points) {
			vec.push_back(p.x);
		}
		return vec;
	}
	vector<float> Yvals() {
		vector<float> vec;
		for (ImVec2 p : this->points) {
			vec.push_back(p.y);
		}
		return vec;
	}
	vector<float> RawXvals() {
		vector<float> vec;
		for (ImVec2 p : this->raw_data) {
			vec.push_back(p.x);
		}
		return vec;
	}
	vector<float> RawYvals() {
		vector<float> vec;
		for (ImVec2 p : this->raw_data) {
			vec.push_back(p.y);
		}
		return vec;
	}
	//function that takes a draw fragment and gets rid of any backwards movement due to precision errors
	void RemoveBackwards() {
		queue<int> back_que;
		queue<int> front_que;
		bool back = false;
		//dont do anything if the function is mixed.
		//TODO: eventually, a warning will pop up along the lines of "loss of precision due to performing *action* on mixed function
		//but since this is only an issue with precision problems, its not a real error, just a warning.
		if (dir == Mixed)
			return;

		//move all negative X values up to be a flat line
		if (dir == Forward) {
			for (int i = 1; i < points.size(); i++) {
				if (points[i].x < points[i - 1].x) {
					points[i].x = points[i - 1].x;
					if (!back) {
						back_que.push(i);
						back = true;
					}
					continue;
				}
				if (back && points[i].x >= points[back_que.front()].x) {
					front_que.push(i);
					back = false;
					continue;
				}
			}
		}
		else {
			for (int i = 1; i < points.size(); i++) {
				if (points[i].x > points[i - 1].x) {
					points[i].x = points[i - 1].x;
					if (!back) {
						back_que.push(i);
						back = true;
					}
					continue;
				}
				if (back && points[i].x <= points[back_que.front()].x) {
					front_que.push(i);
					back = false;
					continue;
				}
			}
		}
		//if the new data point is not between the front and back of the unaffected drawing, then it needs to be removed
		//this prevents awkward vertical line artifacts
		while (!back_que.empty()) {
			int front_index = points.size() - 1;
			int back_index = back_que.front();
			back_que.pop();
			if (!front_que.empty()) {
				front_index = front_que.front();
				front_que.pop();
			}
			float front_val = points[front_index].y;
			float back_val = points[back_index].y;
			for (int i = back_index; i < front_index; i++) {
				points[i].y = back_val + (i - back_index) * GetSlope(ImVec2(back_index, back_val), ImVec2(front_index, front_val));
			}
		}

	}
	void ClearRegression() {
		reg_degree = -1;
		reg_type = None;

	}
	void ResetData() {
		points.assign(raw_data.begin(), raw_data.end());

	}
	//function for calculating manipulation box points around a fragment
	//array is clockwise, starting top left
	array<ImVec2, 8> BoxPoints(bool raw) {
		vector<ImVec2>* data_points = &points;
		if (raw)
			data_points = &raw_data;

		float max_x = (*max_element((*data_points).begin(), (*data_points).end(), [](ImVec2 lhs, ImVec2 rhs) {return lhs.x < rhs.x ? true : false; })).x;
		float max_y = (*max_element((*data_points).begin(), (*data_points).end(), [](ImVec2 lhs, ImVec2 rhs) {return lhs.y < rhs.y ? true : false; })).y;
		float min_x = (*min_element((*data_points).begin(), (*data_points).end(), [](ImVec2 lhs, ImVec2 rhs) {return lhs.x < rhs.x ? true : false; })).x;
		float min_y = (*min_element((*data_points).begin(), (*data_points).end(), [](ImVec2 lhs, ImVec2 rhs) {return lhs.y < rhs.y ? true : false; })).y;
		array<ImVec2, 8> arr;
		arr[0] = ImVec2(min_x, min_y);
		arr[1] = ImVec2((min_x + max_x) / 2, min_y);
		arr[2] = ImVec2(max_x, min_y);
		arr[3] = ImVec2(max_x, (max_y + min_y) / 2);
		arr[4] = ImVec2(max_x, max_y);
		arr[5] = ImVec2((min_x + max_x) / 2, max_y);
		arr[6] = ImVec2(min_x, max_y);
		arr[7] = ImVec2(min_x, (max_y + min_y) / 2);
		return arr;
	}
	//stretches data for manipulation
	//returns a pair of vectors, deltas of the first and last point in the changed vector
	pair<ImVec2, ImVec2> StretchData(ImVec2 mouse_delta, ImVec2 control, ImVec2 dimension, bool raw) {
		pair<ImVec2, ImVec2> edge_deltas;
		vector<ImVec2>* data_points = &points;
		if (raw)
			data_points = &raw_data;
		//calculating first and last for pair (i know i do this calc twice but its negligable
		edge_deltas.first = ImVec2(mouse_delta.x * (1 - abs((*data_points)[0].x - control.x) / dimension.x),
			mouse_delta.y * (1 - abs((*data_points)[0].y - control.y) / dimension.y));
		edge_deltas.second = ImVec2(mouse_delta.x * (1 - abs((*data_points).back().x - control.x) / dimension.x),
			mouse_delta.y * (1 - abs((*data_points).back().y - control.y) / dimension.y));

		//calculating and applying changes to every vector
		float dx, dy;
		for (ImVec2& point : *data_points) {
			dx = mouse_delta.x * (1 - abs(point.x - control.x) / dimension.x);
			dy = mouse_delta.y * (1 - abs(point.y - control.y) / dimension.y);
			point.x += dx;
			point.y += dy;
		}
		return edge_deltas;
	}
};

//fragment equality overload
bool operator ==(const DrawFragment& lhs, const DrawFragment& rhs) {
	return (lhs.points == rhs.points &&
		lhs.raw_data == rhs.raw_data &&
		lhs.dir == rhs.dir &&
		lhs.div == rhs.div &&
		lhs.reg_type == rhs.reg_type &&
		lhs.reg_degree == rhs.reg_degree
		);
}



//state node for the state tree. Currently I have just undo/redo stack, but I may use an undo/redo tree in the future so im keeping this struct
struct StateNode {
	vector<DrawFragment> fragments; //draw data to be saved
	StateNode* parent;
	vector<StateNode*> children;
	StateNode(vector<DrawFragment> fragments, StateNode* parent) {
		this->fragments = fragments;
		this->parent = parent;
	}
};



//checks for mixed directional data
bool ContainsMixed(vector<DrawFragment>& fragments) {
	for (auto f : fragments) {
		if (f.dir & Mixed)
			return true;
	}
	return false;
}

//removes directionality
void SetAllMixed(vector<DrawFragment>& fragments) {
	for (auto& f : fragments)
		f.dir = Mixed;
}

void ResetAll(vector<DrawFragment>& fragments) {
	for (auto& f : fragments)
		f.ResetData();
}

//short function to check if all fragments are empty
bool IsEmptyDrawing(vector<DrawFragment>& fragments) {
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	for (auto& f : fragments)
		if (f.points.size() > 0)
			return false;
	return true;
}

//warning function for too high a regression degree
void RegDegreeError(bool high) {
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	if (ImGui::BeginPopupModal("Reg Degree Error", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
		if (high)
			ImGui::Text("Error: Regression degree too high for floating points to handle");
		else
			ImGui::Text("Error: Regression degree too low");
		if (ImGui::Button("Ok")) {
			ImGui::CloseCurrentPopup();
		}
		ImGui::EndPopup();
	}
}

//error function for an action on an empty drawing
void EmptyCanvasError() {

	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	if (ImGui::BeginPopupModal("Error Empty", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
		ImGui::Text("Error: Cannot perform action on empty canvas");
		if (ImGui::Button("Close"))
			ImGui::CloseCurrentPopup();
		ImGui::EndPopup();
	}
}

void SizeError() {
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
	if (ImGui::BeginPopupModal("Size Error", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
		ImGui::Text("Error: Cannot perform action, fragment size too small");
		if (ImGui::Button("Close"))
			ImGui::CloseCurrentPopup();
		ImGui::EndPopup();
	}
}

//splits a fragment in two, divided at point_index and its corresponding raw_index
//returns the # of splits performed
bool SplitFragments(vector<DrawFragment>& fragments, int frag_index, int point_index, Direction dir1, Direction dir2, DivType divtype) {
	vector<ImVec2> points1, points2, raw1, raw2;
	if (point_index >= fragments[frag_index].points.size() || point_index < 0)
		return false;
	DivType base_div = fragments[frag_index].div;
	int raw_index = floor(static_cast<double>(point_index) / fragments[frag_index].points.size() * fragments[frag_index].raw_data.size());
	//if we try to split at an index past the end of the fragment, we just dont do anything


	for (int i = 0; i <= point_index; i++)
		points1.push_back(fragments[frag_index].points[i]);
	for (int j = point_index; j < fragments[frag_index].points.size(); j++)
		points2.push_back(fragments[frag_index].points[j]);

	for (int i = 0; i <= raw_index; i++)
		raw1.push_back(fragments[frag_index].raw_data[i]);
	for (int j = raw_index; j < fragments[frag_index].raw_data.size(); j++)
		raw2.push_back(fragments[frag_index].raw_data[j]);

	DrawFragment frag1(points1, dir1, divtype);
	frag1.raw_data = raw1;
	DrawFragment frag2(points2, dir2, base_div);
	frag2.raw_data = raw2;
	fragments.erase(fragments.begin() + frag_index);
	fragments.insert(fragments.begin() + frag_index, frag1);
	fragments.insert(fragments.begin() + frag_index + 1, frag2);
	return true;
}

//Helper function to split a fragment multiple times at once
//points and raws need to be the same size, dirs needs to be one size larger
//points need to be sorted in increasing order (otherwise this would be a massive headache)
int SplitMultiple(vector<DrawFragment>& fragments, int frag_index, vector<int> points, vector<Direction> dirs, vector<DivType> divs) {
	int splits = 0;
	int back = 0;
	bool success = false;
	if (points.size() != (dirs.size() - 1))
		return 0;
	for (int i = 0; i < points.size(); i++) {
		success = SplitFragments(fragments, frag_index + splits, points[i] - back, dirs[i], dirs[i + 1], divs[i]);
		if (success) {
			splits++;
			back = points[i];
		}
	}
	return splits;
}



// Merges a fragment into the fragment beyond it in the fragment vector.
// frag_index is the base fragment to be merged into the following fragment
// returns true if the merge was successful
bool MergeFragments(vector<DrawFragment>& fragments, int frag_index) {
	if (frag_index + 1 >= fragments.size() || fragments.size() < 2)
		return false;
	//create references to the data then combine it
	vector<ImVec2>& points1 = fragments[frag_index].points;
	vector<ImVec2>& raw1 = fragments[frag_index].raw_data;
	vector<ImVec2>& points2 = fragments[frag_index + 1].points;
	vector<ImVec2>& raw2 = fragments[frag_index + 1].raw_data;
	points1.insert(points1.end(), points2.begin() + 1, points2.end());
	raw1.insert(raw1.end(), raw2.begin(), raw2.end());
	//swap so we have the combined data in the fragment we keep
	points1.swap(points2);
	raw1.swap(raw2);
	//fixing direction in case they are mixed
	if (fragments[frag_index].dir != fragments[frag_index + 1].dir)
		fragments[frag_index + 1].dir = Mixed;
	//get rid of old fragment
	fragments.erase(fragments.begin() + frag_index);
	//clear the regression data of the fragment because it is no longer regressed
	fragments[frag_index].ClearRegression();

	return true;
}

//helper function for merging multiple fragments together
//returns an integer of the number of successful merges
int MergeMultiple(vector<DrawFragment>& fragments, int begin_index, int end_index, DivType div_condition) {
	int merges = 0;
	int loc = begin_index;
	end_index = min(end_index, static_cast<int>(fragments.size() - 1));
	while (loc <= end_index - merges) {
		//Dont merge if they dont meet the condition
		if (!(div_condition & fragments[loc].div))
			loc++;
		else if (MergeFragments(fragments, loc))
			merges++;
		else
			loc++;
		//TODO finish this
	}
	return merges;
}

//proximity checking function, checking if the cursor is close to a given point
bool WithinReach(ImVec2 mousePos, ImVec2 targetPos, float radius) {
	if (abs(mousePos.x - targetPos.x) < radius && abs(mousePos.y - targetPos.y) < radius)
		return true;
	else
		return false;
}

//simple slope function I use on occasion
float GetSlope(ImVec2 p1, ImVec2 p2) {
	float yDiff = p2.y - p1.y;
	float xDiff = p2.x - p1.x;
	return yDiff / xDiff;
}

//linear search to find the index of a vector val point
//implemented linearly since conditioning/sorting the data would be more expensive than a linear search
int PointIndex(vector<ImVec2> data, ImVec2 point) {
	int count = 0;
	while (count < data.size() && !(data[count] == point))
		count++;
	if (count == data.size())
		return -1;
	return count;
}

//distance between two points
float Distance(ImVec2 p1, ImVec2 p2) {
	return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

//simple linear search for index at which distance from fragment to a point is a minimum
//returns both the index aand the value of the distance as a pair
pair<int, float> MinDistance(DrawFragment& frag, ImVec2 point) {
	vector<ImVec2>& points = frag.points;
	if (points.size() == 0)
		return { 0,0 };
	int min_index = 0;
	float min_val = Distance(points[0], point);
	float temp = 0;
	for (int i = 0; i < points.size(); i++) {
		temp = Distance(points[i], point);
		if (temp < min_val) {
			min_val = temp;
			min_index = i;
		}
	}
	return { min_index, min_val };
}

//linear interpolation for increasing regression accuracy (linearly spaces the data points)
vector<ImVec2> LinearInterpolation(vector<ImVec2>& points, float spacing) {
	vector<ImVec2> interPoints;
	if (points.size() == 0)
		return interPoints;
	ImVec2 p1, p2, curr;
	float distance, slope, theta, remainder;
	remainder = 0;
	interPoints.push_back(points[0]);
	for (int i = 0; i < points.size() - 1; i++) {
		p1 = points[i];
		p2 = points[i + 1];
		distance = Distance(p1, p2);
		slope = GetSlope(p1, p2);
		theta = atan(slope);
		curr = points[i];
		while (distance >= (spacing - remainder)) {
			distance -= (spacing - remainder);
			if (p2.x > p1.x)
				interPoints.push_back(ImVec2(curr.x + ((spacing - remainder) * cos(theta)), curr.y + ((spacing - remainder) * sin(theta))));
			else if (p2.x == p1.x)
				if (p2.y > p1.y)
					interPoints.push_back(ImVec2(curr.x, curr.y + spacing - remainder));
				else
					interPoints.push_back(ImVec2(curr.x, curr.y - spacing - remainder));
			else
				interPoints.push_back(ImVec2(curr.x - ((spacing - remainder) * cos(theta)), curr.y - ((spacing - remainder) * sin(theta))));
			remainder = 0;
			curr = interPoints.back();
		}
		remainder += distance;
	}
	interPoints.push_back(points.back());
	return interPoints;
}

//doesn't actually add to the draw list, just draws it into a vector
vector<ImVec2> DrawBezier(vector<pair<float, float>> cpoints, int resolution) {
	int degree = cpoints.size() - 1;
	vector<ImVec2> draw_points;
	pair<float, float> point;
	for (int i = 0; i <= resolution; i++) {
		point = BezValue(cpoints, i / static_cast<float>(resolution), degree);
		draw_points.push_back(ImVec2(point.first, point.second));
	}
	return draw_points;
}

//doesnt add to draw list, just draws regression into a vector
//function only works on a set of forwards/backwards data, error checking for this is on the side of the main function
vector<ImVec2> DrawPolynomial(vector<float> coeffs, vector<float> xvals, int resolution) {
	float yval = 0; float xval = 0;
	float delta = xvals.back() - xvals[0];
	vector<ImVec2> draw_points;
	for (int i = 0; i <= resolution; i++) {
		yval = 0;
		xval = xvals[0] + delta * i / static_cast<float>(resolution);
		for (int j = 0; j < coeffs.size(); j++) {
			yval += coeffs[j] * pow(xval, j);
		}
		draw_points.push_back(ImVec2(xval, yval));
	}
	return draw_points;
}


//function that finds local mins and maxes
//returns the number of mins/maxes
int DivideMinsMaxes(vector<DrawFragment>& fragments, int precision, std::function<bool(DrawFragment)> lambda) {
	int splits = 0;
	int size = fragments.size();
	for (int f = 0; f < size; f++) {
		DrawFragment copy = fragments[f + splits];
		//ignore min/max if it fits this function criteria
		if (lambda(copy))
			continue;
		copy.RemoveBackwards();
		vector<float> y = copy.Yvals();
		vector<float> x = copy.Xvals();
		vector<int> indexes;
		vector<Direction> directions;
		vector<DivType> div_types;

		Direction dir = fragments[f + splits].dir;
		directions.push_back(dir);
		int j = 0;
		int left = 0;
		int right = 0;
		for (int i = 1; i < copy.points.size() - 1; i++) {
			//check mins
			if (y[i] >= y[i - 1] && y[i] >= y[i + 1]) { //checks if point is potential max
				j = i;
				left = 0;
				right = 0;
				while (y[j - left] <= y[j] && y[j + right] <= y[j]) {
					if (y[j] - y[j - left] >= precision / 2.0 && (fabs(x[j] - x[j - left]) >= precision / 2.0 || (j - left) == 0)) //threshold criteria left
						if (y[j] - y[j + right] >= precision / 2.0 && (fabs(x[j] - x[j + right]) >= precision / 2.0 || (j + right) == y.size() - 1)) { //threshold criteria right
							//bonus loop that helps find the middle *visual* min
							int mid_left = 0;
							int mid_right = 0;
							for (mid_right = 0; mid_right < right; mid_right++) {
								if (abs(y[j] - y[j + mid_right]) > .25)
									break;
							}
							for (mid_left = 0; mid_left < left; mid_left++) {
								if (abs(y[j] - y[j - mid_left]) > .25)
									break;
							}
							j += (mid_left + mid_right) / 2;
							indexes.push_back(PointIndex(fragments[f + splits].points, ImVec2(x[j], y[j])));
							directions.push_back(dir);
							div_types.push_back(Extrema);
							i += right;
							break;

						}
					if (j - left > 0)
						left++;
					if (j + right < y.size() - 1)
						right++;
					if (j - left <= 0 && j + right >= y.size() - 1)
						break;
				}
			}
		}
		for (int i = 1; i < copy.points.size() - 1; i++) {
			//check Maxes
			if (y[i] <= y[i - 1] && y[i] <= y[i + 1]) { //checks if point is potential max
				j = i;
				left = 0;
				right = 0;
				while (y[j - left] >= y[j] && y[j + right] >= y[j]) {
					if (y[j] - y[j - left] <= -precision / 2.0 && (fabs(x[j] - x[j - left]) >= precision / 2.0 || (j - left) == 0)) //threshold criteria left
						if (y[j] - y[j + right] <= -precision / 2.0 && (fabs(x[j] - x[j + right]) >= precision / 2.0 || (j + right) == y.size() - 1)) { //threshold criteria right
							//bonus loop that helps find the middle *visual* min
							int mid_left = 0;
							int mid_right = 0;
							for (mid_right = 0; mid_right < right; mid_right++) {
								if (abs(y[j] - y[j + mid_right]) > .25)
									break;
							}
							for (mid_left = 0; mid_left < left; mid_left++) {
								if (abs(y[j] - y[j - mid_left]) > .25)
									break;
							}
							j += (mid_left + mid_right) / 2;
							indexes.push_back(PointIndex(fragments[f + splits].points, ImVec2(x[j], y[j])));
							directions.push_back(dir);
							div_types.push_back(Extrema);
							i += right;
							break;

						}
					if (j - left > 0)
						left++;
					if (j + right < y.size() - 1)
						right++;
					if (j - left <= 0 && j + right >= y.size() - 1)
						break;
				}
			}
		}
		std::sort(indexes.begin(), indexes.end());
		splits += SplitMultiple(fragments, f + splits, indexes, directions, div_types);
	}
	return splits;
}





void undo(stack<StateNode*>& redo_stack, StateNode*& base_node, vector<DrawFragment>& fragments, vector<ImVec2>& points, bool& proximity_flag) {
	if (base_node->parent == nullptr)
		return;
	redo_stack.push(base_node);
	base_node = base_node->parent;
	fragments = base_node->fragments;
	points = fragments.back().points;

	//if this is after an erase, allow for new drawing anywhere (if not intensive, maybe just set this flag to true whenever size is zero)
	if (points.size() == 0)
		proximity_flag = true;
}

void redo(stack<StateNode*>& redo_stack, StateNode*& base_node, vector<DrawFragment>& fragments, vector<ImVec2>& points, bool& proximity_flag) {
	if (!redo_stack.empty()) {
		base_node = redo_stack.top();
		redo_stack.pop();
		fragments = base_node->fragments;
		points = fragments.back().points;

		if (points.size() == 0)
			proximity_flag = true;

	}
}

//re-usable widgets
void Precision(int& precision) {
	ImGui::PushItemWidth(100);
	ImGui::InputInt("input precision", &precision);
	ImGui::PopItemWidth();
}
void Depth(int& depth) {
	ImGui::PushItemWidth(100);
	ImGui::InputInt("input depth", &depth);
	ImGui::PopItemWidth();
}
void Degree(int& degree) {
	ImGui::PushItemWidth(100);
	ImGui::InputInt("input degree", &degree);
	ImGui::PopItemWidth();
}
//TODO add "don't show me this again" checkbox and add it to every warning (but not error obv)

int main()
{
	//Setting up GLFW connections
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	GLFWwindow* window = glfwCreateWindow(1280, 720, "Wave Draw", NULL, NULL);
	if (window == NULL) {
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	//IMGUI prep stuff
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 430"); //using the shaders for the latest supported openGL version 4.3
	bool show_demo_window = true;
	bool show_another_window = true;
	bool waveDrawOpen = true;
	bool startup = true;
	ImVec4 clear_color = ImVec4(.45f, .55f, .60f, 1.00f);

	//Main loop
	while (!glfwWindowShouldClose(window)) {
		//poll inputs then update to a new frame
		glfwPollEvents();
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		//DEMO WINDOW
		if (show_demo_window)
			ImGui::ShowDemoWindow(&show_demo_window);

		//MY WINDOW
		//variable and structure initialization (everything here is accessible by every GUI window/menu/etc. scope)
		WD_Window current = Drawing;
		static vector<ImVec2> points;
		static vector<DrawFragment> fragments;
		static StateNode* base_node = new StateNode(fragments, nullptr); //modifications apply to THIS node
		static stack<StateNode*> redo_stack;
		static bool proximity_flag = true; //is cursor hovering near drawing point?

		//anything to be done before the GUI actually starts
		if (startup) {
			fragments.push_back(DrawFragment(points, Mixed, Other));
			base_node->fragments = fragments;
			startup = false;

		}

		//TODO main menu
		if (ImGui::BeginMainMenuBar()) {
			if (ImGui::BeginMenu("File")) {
				if (current == Drawing) {
					if (ImGui::MenuItem("Save", "CTRL-S")) {}
					if (ImGui::MenuItem("Save-As Wave Drawing")) {}
				}
				else if (current == Shaping || current == Superposition) {
					if (ImGui::MenuItem("Save", "CTRL-S")) {}
					if (ImGui::MenuItem("Save-As Wave Vector")) {}
				}
				else if (current == Playback) {
					if (ImGui::MenuItem("Save", "CTRL-S")) {}
					if (ImGui::MenuItem("Save-As Audio File")) {}
				}
				else { //not sure what to do here
					if (ImGui::MenuItem("Save", "CTRL-S")) {}
					if (ImGui::MenuItem("Save-As")) {}
				}
				ImGui::EndMenu();
			}
			if (ImGui::BeginMenu("Edit")) {
				if (ImGui::MenuItem("Undo", "CTRL-Z")) {
					undo(redo_stack, base_node, fragments, points, proximity_flag);
				}
				if (ImGui::MenuItem("Redo", "CTRL-Y")) {
					redo(redo_stack, base_node, fragments, points, proximity_flag);
				}
				ImGui::EndMenu();
			}

			ImGui::EndMainMenuBar();
		}

		//if(current == Drawing) use this once multiple windows are being used
		if (waveDrawOpen)
		{
			//Drawing window variable & lambda initialization
			static ImVec2 origin = ImVec2(0, 0); //initialize origin as static var so we can pass it between frames and calc mouse position at beginning of frame								 	
			static ImVec2 canvas_p0 = ImVec2(0, 0);      // ImDrawList API uses screen coordinates!;
			static ImVec2 canvas_sz = ImVec2(0, 0);   // Resize canvas to what's available
			static ImVec2 canvas_p1 = ImVec2(0, 0);
			const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);
			const ImVec2 old_mouse_pos(io.MousePos.x - origin.x - io.MouseDelta.x, io.MousePos.y - origin.y - io.MouseDelta.y);
			static bool is_hovered = false; // Hovered
			static bool is_active = false;   // Held
			static bool opt_enable_grid = true;
			static bool rc_popup = false;
			static bool reshaping = false;
			static bool manip = false;
			static bool disable = false;
			static bool direction_colors = true;
			static bool drawing = false; //is user currently drawing?
			static bool new_state = false; // flag for creating a new state node
			static ImVec2 scrolling(0.0f, 0.0f);
			static int precision = 5;
			static int state_count = 0;
			static int clicked_frag = -1;
			static int clicked_control = -1;
			static int rc_frag = 0;




			//for now I need to include arbitrary param, later may use function templates and check to make sure the lambda is a predicate
			std::function<bool(DrawFragment)> check_mixed = [](DrawFragment frag) { return frag.dir == Mixed ? true : false; };
			std::function<bool(DrawFragment)> return_true = [](DrawFragment frag) { return true; };
			std::function<bool(DrawFragment)> return_false = [](DrawFragment frag) { return false; };

			//points vector is only what can be drawn on, so this is the end fragment of the drawing
			fragments.back().points = points;

			//check if a new state needs to be saved
			if (new_state) {
				//if the data doesn't change, we don't create a new state
				if (!(fragments == base_node->fragments)) {
					//create new node here
					StateNode* new_node = new StateNode(fragments, base_node);
					//push back the new node
					base_node->children.push_back(new_node);
					base_node = new_node;
					//get rid of redo stack
					while (!redo_stack.empty()) {
						redo_stack.pop();
					}
					state_count++;
				}
				new_state = false;
			}

			//setting up main window
			ImGuiWindowFlags waveDrawFlags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings;
			//waveDrawFlags |= ImGuiWindowFlags_NoTitleBar  | ImGuiWindowFlags_NoCollapse; uncomment when done with demo file
			const ImGuiViewport* viewport = ImGui::GetMainViewport();
			ImGui::SetNextWindowPos(viewport->WorkPos);
			ImGui::SetNextWindowSize(viewport->WorkSize);

			//hotkey functions
			auto undo_lambda = []() {undo(redo_stack, base_node, fragments, points, proximity_flag); };
			auto redo_lambda = []() {redo(redo_stack, base_node, fragments, points, proximity_flag); };
			//TODO add save hotkeys (not sure how to handle saving via multiple options, but maybe the lambda function can include a switch statement
			//     that calls different save/save as functions depending on what the current window is (current = Drawing, current = Modulation, etc.)
			static bool no_keys = false;
			static int key_count = 0;
			map<ImGuiKey, function<void()>> hotkeys = { {ImGuiKey_Z, undo_lambda}, {ImGuiKey_Y, redo_lambda} };
			map<ImGuiKey, function<void()>>::iterator it;
			if (io.KeyCtrl) {
				key_count = 0;
				for (it = hotkeys.begin(); it != hotkeys.end(); it++) {
					if (ImGui::IsKeyDown(it->first)) {
						if (no_keys)
							it->second();
						key_count++;
					}
				}
				if (key_count == 0)
					no_keys = true;
				else
					no_keys = false;
			}
			else
				no_keys = false;


			//Start the window frame
			ImGui::Begin("WaveDraw", NULL, waveDrawFlags);
			//creating draw list and splitting it into render sections
			ImDrawList* draw_list = ImGui::GetWindowDrawList();
			ImDrawListSplitter draw_splitter;
			draw_splitter.Split(draw_list, 3);

			// I want all the draw list calls before my clipping rect to be rendered after the clipping rect, so I make it like this
			draw_splitter.SetCurrentChannel(draw_list, 1);

			//modification tab bar
			if (ImGui::BeginTabBar("Wave Draw Tabs"))
			{
				if (ImGui::BeginTabItem("Drawing Station")) {
					reshaping = false;
					if (ImGui::Button("Erase Drawing")) {
						if (points.size() != 0 || fragments.size() != 1)
							new_state = true;
						points.clear();
						fragments.back().raw_data.clear();
						fragments.clear();
						fragments.push_back(DrawFragment(points, Mixed, Other));
						proximity_flag = true;
						disable = false;
					}
					ImGui::EndTabItem();
				}
				if (ImGui::BeginTabItem("Breakpoints")) {
					reshaping = false;
					//TODO merge all directionality fragments before recalculating directionality (also throw up a warning about this with a dont show again box)
					if (ImGui::Button("Calculate Directionality")) {
						//divides drawing up into forwards and backwards function
						int fsize = fragments.size();
						int total_broken = 0;
						for (int f = 0; f < fsize + total_broken; f++) {
							vector<ImVec2> frag_points = fragments[f].points;
							int size = frag_points.size();
							int back = 0;
							int erased = 0;
							int broken = 0;
							bool window = false;
							bool inv_window = false;
							for (int i = 1; i < size; i++) {
								if (window) {
									if (frag_points[back].x - frag_points[i].x > precision) {
										//a little hack that catches fragments that only go backwards
										if (back == 0) {
											fragments[f].dir = Backward;
										}
										else {
											SplitFragments(fragments, f + broken, back - erased, Forward, Backward, ChangeDir);
											broken++;
										}
										window = false;
										inv_window = true;
										erased = back;
										back = i;

									}
									else if (frag_points[i].x >= frag_points[back].x) {
										back = i;
									}
								}
								else if (inv_window) {
									if (frag_points[i].x - frag_points[back].x > precision) {
										SplitFragments(fragments, f + broken, back - erased, Backward, Forward, ChangeDir);
										inv_window = false;
										erased = back;
										back = i;
										broken++;
									}
									else if (frag_points[i].x <= frag_points[back].x) {
										back = i;
									}
								}
								else if (frag_points[i].x <= frag_points[i - 1].x) {
									back = i - 1;
									window = true;
								}

							}
							//a little hack that catches fragments that dont change direction at all
							if (broken == 0 && frag_points.size() > 0) {
								if (frag_points.back().x >= frag_points[0].x)
									fragments[f].dir = Forward;
								else
									fragments[f].dir = Backward;
							}
							total_broken += broken;
						}
						new_state = true;
						points = fragments.back().points;
					} ImGui::SameLine();
					if (ImGui::Button("Clear Directionality")) {
						//merge changed directions
						MergeMultiple(fragments, 0, fragments.size() - 1, (ChangeDir));
						SetAllMixed(fragments);
						new_state = true;
						points = fragments.back().points;
					} ImGui::SameLine();

					//direction color checkbox
					ImGui::Checkbox("Direction Colors", &direction_colors); ImGui::SameLine();

					//divides the function up by direction (will eventually open up menus to do regressions, and further function division)
					if (ImGui::Button("Find Extrema")) {
						if (IsEmptyDrawing(fragments)) {
							if (!(ImGui::IsPopupOpen("Error Empty", ImGuiPopupFlags_AnyPopupId)))
								ImGui::OpenPopup("Error Empty", ImGuiWindowFlags_AlwaysAutoResize);
						}
						else if (ContainsMixed(fragments)) {
							if (!(ImGui::IsPopupOpen("Error", ImGuiPopupFlags_AnyPopupId)))
								ImGui::OpenPopup("MixedExtrema");
						}
						else {
							DivideMinsMaxes(fragments, precision, return_false);
							new_state = true;
							points = fragments.back().points;
						}

					} ImGui::SameLine();



					//mixed extrema popup
					ImVec2 center = ImGui::GetMainViewport()->GetCenter();
					ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
					if (ImGui::BeginPopupModal("MixedExtrema", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
						ImGui::Text("Mixed fragments found, do you want to skip Min/Max calculation for these?");
						ImGui::Separator();
						if (ImGui::Button("Skip")) {
							DivideMinsMaxes(fragments, precision, check_mixed);
							new_state = true;
							points = fragments.back().points;
							ImGui::CloseCurrentPopup();
						}
						ImGui::SameLine();
						if (ImGui::Button("Dont Skip")) {
							DivideMinsMaxes(fragments, precision, return_false);
							new_state = true;
							points = fragments.back().points;
							ImGui::CloseCurrentPopup();
						}
						ImGui::EndPopup();
					}
					//empty popup
					EmptyCanvasError();

					if (ImGui::Button("Clear Extrema")) {
						MergeMultiple(fragments, 0, fragments.size() - 1, Extrema);
						new_state = true;
						points = fragments.back().points;
					} ImGui::SameLine();

					//merges all functions in case the user wants to reset the function back to normal
					if (ImGui::Button("Remove All")) {
						MergeMultiple(fragments, 0, fragments.size() - 1, (Manual | Extrema | ChangeDir | Other));
						new_state = true;
						points = fragments.back().points;
					} ImGui::SameLine();

					//Precision Input
					Precision(precision);

					//breakpoint sliding (may make this into a function)
					static bool bp_clicked = false;
					static bool bp_dragging = false;
					static int bp_index = 0;
					static int temp_index = 0;
					static float bp_dist = FLT_MAX;
					static float temp_dist = FLT_MAX;
					static int bp_frag = 0;
					if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && points.size() != 0) {
						for (int i = 0; i < fragments.size() - 1; i++) {
							if (WithinReach(mouse_pos_in_canvas, fragments[i].points.back(), 5)) {
								bp_clicked = true;
								bp_dragging = true;
								MergeFragments(fragments, i);
								points = fragments.back().points;
								break;
							}
						}
					}
					if (bp_clicked) {
						if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
							//TODO write a better algorithm for this when you have the time
							// (maybe some type of binary search, but has to be linear at some point idk)
							bp_index = 0;
							temp_index = 0;
							bp_dist = FLT_MAX;
							temp_dist = FLT_MAX;
							for (int i = 0; i < fragments.size(); i++) {
								tie(temp_index, temp_dist) = MinDistance(fragments[i], mouse_pos_in_canvas);
								if (temp_dist < bp_dist) {
									bp_dist = temp_dist;
									bp_index = temp_index;
									bp_frag = i;
								}
							}
							ImVec2 bp_point = fragments[bp_frag].points[bp_index];
							draw_list->AddCircle(ImVec2(origin.x + bp_point.x, origin.y + bp_point.y), 8, IM_COL32(0, 255, 155, 255), 20, 2.0);
						}
						else {
							SplitFragments(fragments, bp_frag, bp_index, fragments[bp_frag].dir, fragments[bp_frag].dir, Other);
							points = fragments.back().points;
							bp_clicked = false;
							bp_dragging = false;
							new_state = true;

						}
					}

					//click to add breakpoint
					if (!bp_dragging && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left) && points.size() != 0) {
						float dbl_dist = FLT_MAX;
						float temp_dist_dbl = FLT_MAX;
						int temp_index_dbl = 0;
						int dbl_index = 0;
						int dbl_frag = 0;
						for (int i = 0; i < fragments.size(); i++) {
							//find the closest spot
							tie(temp_index_dbl, temp_dist_dbl) = MinDistance(fragments[i], mouse_pos_in_canvas);
							if (temp_dist_dbl < dbl_dist) {
								dbl_dist = temp_dist_dbl;
								dbl_index = temp_index_dbl;
								dbl_frag = i;
							}
						}
						//check if its close enough to be a double click and split it
						vector<ImVec2>& dbl_points = fragments[dbl_frag].points;
						if (dbl_dist < 5) {
							SplitFragments(fragments, dbl_frag, dbl_index, fragments[dbl_frag].dir, fragments[dbl_frag].dir, Manual);
							points = fragments.back().points;
							new_state = true;
						}
					}
					if (!bp_dragging && ImGui::IsMouseClicked(ImGuiMouseButton_Right) && points.size() != 0) {
						float rc_dist = FLT_MAX;
						float temp_dist_rc = FLT_MAX;
						rc_frag = 0;
						for (int i = 0; i < fragments.size(); i++) {
							//find the closest spot
							temp_dist_rc = Distance(fragments[i].points.back(), mouse_pos_in_canvas);
							if (temp_dist_rc < rc_dist) {
								rc_dist = temp_dist_rc;
								rc_frag = i;
							}
						}
						//if its close to a breakpoint then we will open a context menu for deleting said breakpoint
						if (rc_dist < 5) {
							rc_popup = true;
						}
					}
					//TODO add delete key and/or backspace key hotkey that deletes a selected breakpoint


					ImGui::EndTabItem();
				}
				if (ImGui::BeginTabItem("Reshaping")) {
					static array<ImVec2, 8> box_points;
					static array<ImVec2, 8> box_points_raw;


					if (clicked_frag != -1 && points.size() != 0) {
						box_points = fragments[clicked_frag].BoxPoints(false); //array of control points
						box_points_raw = fragments[clicked_frag].BoxPoints(true); //array of control points
						for (int i = 0; i < box_points.size(); i++) {
							ImVec2 cp = box_points[i];
							draw_list->AddCircle(ImVec2(origin.x + cp.x, origin.y + cp.y), 5, IM_COL32(255, 255, 255, 255), 20, 2.0);
						}
						//TODO draw dashes here
					}
					ImGui::Checkbox("Frag_Select", &reshaping); ImGui::SameLine();
					//logic for selecting a fragment or manipulation point
					if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && reshaping && is_hovered && points.size() != 0) {
						int clickcount = 0;
						if (clicked_frag != -1) {
							int controlcount = 0;
							for (int i = 0; i < box_points.size(); i++) {
								if (WithinReach(mouse_pos_in_canvas, box_points[i], 7)) {
									clicked_control = i;
									controlcount++;
									clickcount++;
									break;
								}
							}
							if (controlcount == 0)
								clicked_control = -1; //not sure if necessary, mainly a failsafe in case clicked_control isn't -1 by default
						}
						if (clicked_control == -1) {
							//checking click fragments
							for (int i = 0; i < fragments.size(); i++) {
								pair<int, float> dist_data = MinDistance(fragments[i], mouse_pos_in_canvas);
								//things to do once a fragment is clicked
								if (dist_data.second < 10) {
									clicked_frag = i;
									clickcount++;
									break;
								}
							}
						}

						if (clickcount == 0) {
							clicked_frag = -1;
							clicked_control = -1;
						}
					}
					//manipulation switch statement determines how to stretch/squish screen
					//TODO add limits to how smooshed the screen can be
					//if height or width is <10 pixels, mouse delta is set to zero if its in smooshing direction
					if (clicked_control != -1 && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
						ImVec2 dimensions = ImVec2(box_points[2].x - box_points[0].x, box_points[4].y - box_points[2].y);
						ImVec2 d_first, d_last, d_first_raw, d_last_raw, translate;
						translate = io.MouseDelta;
						bool minx, miny;
						switch (clicked_control) {
						case 0: //top left
							translate = mouse_pos_in_canvas - box_points[0];
							minx = (dimensions.x - 5) <= abs(translate.x);
							miny = (dimensions.y - 5) <= abs(translate.y);
							if (minx && (translate.x > 0 || mouse_pos_in_canvas.x > box_points[0].x)) {
								translate.x = dimensions.x - 5;
							}
							if (miny && (translate.y > 0 || mouse_pos_in_canvas.y > box_points[0].y)) {
								translate.y = dimensions.y - 5;
							}
								tie(d_first, d_last) = fragments[clicked_frag].StretchData(translate, box_points[0], dimensions, false);
								tie(d_first_raw, d_last_raw) = fragments[clicked_frag].StretchData(translate, box_points_raw[0], dimensions, true);
							break;
						case 1: //top
							translate = mouse_pos_in_canvas - box_points[1];
							minx = (dimensions.x - 5) <= abs(translate.x);
							miny = (dimensions.y - 5) <= abs(translate.y);
							translate.x = 0;
							if (miny && (translate.y > 0 || mouse_pos_in_canvas.y > box_points[1].y))
								translate.y = dimensions.y - 5;
							tie(d_first, d_last) = fragments[clicked_frag].StretchData(translate, box_points[1], dimensions, false);
							tie(d_first_raw, d_last_raw) = fragments[clicked_frag].StretchData(translate, box_points_raw[1], dimensions, true);
							break;
						case 2: // top right
							translate = mouse_pos_in_canvas - box_points[2];
							minx = (dimensions.x - 5) <= abs(translate.x);
							miny = (dimensions.y - 5) <= abs(translate.y);
							if (minx && (translate.x < 0 || mouse_pos_in_canvas.x < box_points[2].x))
								translate.x = -(dimensions.x - 5);
							if (miny && (translate.y > 0 || mouse_pos_in_canvas.y > box_points[2].y))
								translate.y = (dimensions.y - 5);
							tie(d_first, d_last) = fragments[clicked_frag].StretchData(translate, box_points[2], dimensions, false);
							tie(d_first_raw, d_last_raw) = fragments[clicked_frag].StretchData(translate, box_points_raw[2], dimensions, true);
							break;
						case 3: // right
							translate = mouse_pos_in_canvas - box_points[3];
							minx = (dimensions.x - 5) <= abs(translate.x);
							miny = (dimensions.y - 5) <= abs(translate.y);
							if (minx && (translate.x < 0|| mouse_pos_in_canvas.x < box_points[3].x))
								translate.x = -(dimensions.x - 5);
							translate.y = 0;
							tie(d_first, d_last) = fragments[clicked_frag].StretchData(translate, box_points[3], dimensions, false);
							tie(d_first_raw, d_last_raw) = fragments[clicked_frag].StretchData(translate, box_points_raw[3], dimensions, true);
							break;
						case 4: // bottom right
							translate = mouse_pos_in_canvas - box_points[4];
							minx = (dimensions.x - 5) <= abs(translate.x);
							miny = (dimensions.y - 5) <= abs(translate.y);
							if (minx && (translate.x < 0 || mouse_pos_in_canvas.x < box_points[4].x))
								translate.x = -(dimensions.x - 5);
							if (miny && (translate.y < 0 || mouse_pos_in_canvas.y < box_points[4].y))
								translate.y = -(dimensions.y - 5);
							tie(d_first, d_last) = fragments[clicked_frag].StretchData(translate, box_points[4], dimensions, false);
							tie(d_first_raw, d_last_raw) = fragments[clicked_frag].StretchData(translate, box_points_raw[4], dimensions, true);
							break;
						case 5: // bottom
							translate = mouse_pos_in_canvas - box_points[5];
							minx = (dimensions.x - 5) <= abs(translate.x);
							miny = (dimensions.y - 5) <= abs(translate.y);
							translate.x = 0;
							if (miny && (translate.y < 0 || mouse_pos_in_canvas.y < box_points[5].y))
								translate.y = -(dimensions.y - 5);
							tie(d_first, d_last) = fragments[clicked_frag].StretchData(translate, box_points[5], dimensions, false);
							tie(d_first_raw, d_last_raw) = fragments[clicked_frag].StretchData(translate, box_points_raw[5], dimensions, true);
							break;
						case 6: // bottom left
							translate = mouse_pos_in_canvas - box_points[6];
							minx = (dimensions.x - 5) <= abs(translate.x);
							miny = (dimensions.y - 5) <= abs(translate.y);
							if (minx && (translate.x > 0 || mouse_pos_in_canvas.x > box_points[6].x))
								translate.x = dimensions.x - 5;
							if (miny && (translate.y < 0 || mouse_pos_in_canvas.y < box_points[6].y))
								translate.y = -(dimensions.y - 5);
							tie(d_first, d_last) = fragments[clicked_frag].StretchData(translate, box_points[6], dimensions, false);
							tie(d_first_raw, d_last_raw) = fragments[clicked_frag].StretchData(translate, box_points_raw[6], dimensions, true);
							break;
						case 7: // left
							translate = mouse_pos_in_canvas - box_points[7];
							minx = (dimensions.x - 5) <= abs(translate.x);
							miny = (dimensions.y - 5) <= abs(translate.y);
							if (minx && (translate.x > 0 || mouse_pos_in_canvas.x > box_points[7].x))
								translate.x = dimensions.x - 5;
							translate.y = 0;
							tie(d_first, d_last) = fragments[clicked_frag].StretchData(translate, box_points[7], dimensions, false);
							tie(d_first_raw, d_last_raw) = fragments[clicked_frag].StretchData(translate, box_points_raw[7], dimensions, true);
							break;
						}
						//shifting over all the other draw fragments TODO: find a better way to do this 
						for (int i = 0; i < fragments.size(); i++) {
							if (i < clicked_frag) {
								for (ImVec2& p : fragments[i].points)
									p = p + d_first;
								for (ImVec2& p : fragments[i].raw_data)
									p = p + d_first_raw;
							}
							else if (i > clicked_frag) {
								for (ImVec2& p : fragments[i].points)
									p = p + d_last;
								for (ImVec2& p : fragments[i].raw_data)
									p = p + d_last_raw;
							}
						}
						//updating transient points
						points = fragments.back().points;
					}
					else if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
						//TODO STRETCH RAW DATA HERE
						clicked_control = -1;
						points = fragments.back().points; //prevents issues if there is only 1 draw fragment
						new_state = true;
					}


					//regression selection combo
					const char* regressions[] = { "No Regression", "Polynomial", "Bezier" };
					static int item_current_idx = 0;
					ImGui::PushItemWidth(200);
					if (ImGui::BeginCombo("Reg Type", regressions[item_current_idx], None))
					{
						for (int n = 0; n < sizeof(regressions) / sizeof(regressions[0]); n++)
						{
							const bool is_selected = (item_current_idx == n);
							if (ImGui::Selectable(regressions[n], true))
								item_current_idx = n;

							// Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
							if (is_selected)
								ImGui::SetItemDefaultFocus();
						}
						ImGui::EndCombo();
					} ImGui::SameLine();
					ImGui::PopItemWidth();
					//regression activation
					//TODO add error checking for fragment regressions that are too small (must be 1 greater than the degree at least)
					static bool high = true;
					if (regressions[item_current_idx] == "No Regression") {
					}
					else if (regressions[item_current_idx] == "Bezier") {
						static int bez_degree = 3;
						static int bez_depth = 100;
						Degree(bez_degree); ImGui::SameLine();
						Depth(bez_depth); ImGui::SameLine();
						if (ImGui::Button("Regress") && clicked_frag != -1) {
							if (bez_degree > 10) {
								if (!(ImGui::IsPopupOpen("Error Empty", ImGuiPopupFlags_AnyPopupId)))
									ImGui::OpenPopup("Reg Degree Error", ImGuiWindowFlags_AlwaysAutoResize);
								high = true;
							}
							else if (bez_degree < 1) {
								if (!(ImGui::IsPopupOpen("Error Empty", ImGuiPopupFlags_AnyPopupId)))
									ImGui::OpenPopup("Reg Degree Error", ImGuiWindowFlags_AlwaysAutoResize);
								high = false;
							}
							else if (fragments[clicked_frag].raw_data.size() <= bez_degree) {
								if (!(ImGui::IsPopupOpen("Error Empty", ImGuiPopupFlags_AnyPopupId)))
									ImGui::OpenPopup("Size Error");
							}
							else {
								//using a wider spacing because too many data points can slow down the regression
								DrawFragment tempfrag(LinearInterpolation(fragments[clicked_frag].raw_data, 10));
								//through empirical testing, anything above 20 data points usually produces great results
								if (tempfrag.points.size() < 20)
									tempfrag.points = LinearInterpolation(fragments[clicked_frag].raw_data, 10.0 * tempfrag.points.size() / 21);
								fragments[clicked_frag].points = DrawBezier(PolyBezierRegression(tempfrag.Xvals(), tempfrag.Yvals(), bez_depth, bez_degree, true, false), 100);
								fragments[clicked_frag].reg_type = bez_degree == 3 ? BezCubic : bez_degree == 2 ? BezQuadratic : bez_degree == 1 ? BezLinear : BezPolyFit;
								fragments[clicked_frag].reg_degree = bez_degree;
								points = fragments.back().points; //prevents issues if there is only 1 draw fragment
								new_state = true;

							}
						} ImGui::SameLine();
					}
					else if (regressions[item_current_idx] == "Polynomial") {
						static int poly_degree = 1;
						Degree(poly_degree); ImGui::SameLine();
						if (ImGui::Button("Regress") && clicked_frag != -1) {
							if (fragments[clicked_frag].dir & Mixed) {
								if (!(ImGui::IsPopupOpen("Error Empty", ImGuiPopupFlags_AnyPopupId)))
									ImGui::OpenPopup("Mixed Error", ImGuiWindowFlags_AlwaysAutoResize);
							}
							else if (poly_degree > 6) {
								if (!(ImGui::IsPopupOpen("Error Empty", ImGuiPopupFlags_AnyPopupId)))
									ImGui::OpenPopup("Reg Degree Error", ImGuiWindowFlags_AlwaysAutoResize);
								high = true;
							}
							else if (poly_degree < 1) {
								if (!(ImGui::IsPopupOpen("Error Empty", ImGuiPopupFlags_AnyPopupId)))
									ImGui::OpenPopup("Reg Degree Error", ImGuiWindowFlags_AlwaysAutoResize);
								high = false;
							}
							else if (fragments[clicked_frag].raw_data.size() <= poly_degree) {
								if (!(ImGui::IsPopupOpen("Error Empty", ImGuiPopupFlags_AnyPopupId)))
									ImGui::OpenPopup("Size Error");
							}
							else {
								fragments[clicked_frag].points = DrawPolynomial(PolyRegressionFixed(fragments[clicked_frag].RawXvals(), fragments[clicked_frag].RawYvals(), poly_degree), fragments[clicked_frag].RawXvals(), 100);
								fragments[clicked_frag].reg_type = poly_degree == 2 ? Quadratic : poly_degree == 1 ? Linear : PolyFit;
								fragments[clicked_frag].reg_degree = poly_degree;
								points = fragments.back().points;
								new_state = true;
							}

						} ImGui::SameLine();
					}
					//error popups for regressions
					RegDegreeError(high);

					//error popups for data that isnt large enough
					SizeError();

					if (ImGui::BeginPopupModal("Mixed Error", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
						ImGui::Text("Error, cannot perform action on Mixed directional data");
						if (ImGui::Button("OK"))
							ImGui::CloseCurrentPopup();
						ImGui::EndPopup();
					}

					if (ImGui::Button("Reset Curr") && clicked_frag != -1) {
						fragments[clicked_frag].points = fragments[clicked_frag].raw_data;
						fragments[clicked_frag].reg_type = None;
						fragments[clicked_frag].reg_degree = -1;
						points = fragments.back().points;
						new_state = true;
					} ImGui::SameLine();

					if (ImGui::Button("Reset All")) {
						if (!(ImGui::IsPopupOpen("Error Empty", ImGuiPopupFlags_AnyPopupId)))
							ImGui::OpenPopup("Reset Warning");
					}

					ImVec2 center = ImGui::GetMainViewport()->GetCenter();
					ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
					if (ImGui::BeginPopupModal("Reset Warning", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
						ImGui::Text("Warning: Are you sure you want to reset your data? This will reset all regressions");
						if (ImGui::Button("Ok")) {
							ResetAll(fragments);
							new_state = true;
							points = fragments.back().points;
							ImGui::CloseCurrentPopup();
						} ImGui::SameLine();
						if (ImGui::Button("Cancel")) {
							ImGui::CloseCurrentPopup();
						}
						ImGui::EndPopup();
					}
					ImGui::EndTabItem();
				}
				ImGui::EndTabBar();
			}

			draw_splitter.SetCurrentChannel(draw_list, 0);
			draw_list = ImGui::GetWindowDrawList();
			canvas_p0 = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!;
			canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
			if (canvas_sz.x < 50.0f) canvas_sz.x = 50.0f;
			if (canvas_sz.y < 50.0f) canvas_sz.y = 50.0f;
			canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);
			// Draw border and background color
			ImGuiIO& io = ImGui::GetIO();
			draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
			draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));
			origin = ImVec2(canvas_p0.x + scrolling.x, canvas_p0.y + scrolling.y);

			if (!reshaping) {
				//TODO reset clicked_control to -1;
				clicked_frag = -1;
			}
			// This will catch our interactions
			ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
			is_hovered = ImGui::IsItemHovered(); // Hovered
			is_active = ImGui::IsItemActive();   // Held
			//delete breakpoint popup (had to leave it here to fix bugs)
			if (rc_popup) {
				ImGui::OpenPopup("delete_bp_popup", ImGuiWindowFlags_AlwaysAutoResize);
				rc_popup = false;
			}
			if (ImGui::BeginPopup("delete_bp_popup")) {
				if (ImGui::Selectable("Delete")) {
					MergeFragments(fragments, rc_frag);
					ImGui::CloseCurrentPopup();
				}
				if (ImGui::Selectable("Cancel")) {
					ImGui::CloseCurrentPopup();
				}
				ImGui::EndPopup();
			}
			//mouse pos
			ImGui::Text("Mouse X: %f      Mouse Y: %f", mouse_pos_in_canvas.x, mouse_pos_in_canvas.y);
			//state count
			ImGui::Text("State count is: %d", state_count);

			// Check for proximity
			if (points.size() > 0 && !ImGui::IsMouseDown(ImGuiMouseButton_Left) && !reshaping)
			{
				if (WithinReach(mouse_pos_in_canvas, points.back(), 5)) {
					proximity_flag = true;
					glfwSetCursorPos(window, points.back().x + origin.x, points.back().y + origin.y);
				}
				else
					proximity_flag = false;
			}

			//add movement data
			//TODO make sure you cant draw and move a breakpoint at the same time (add another flag)
			if (ImGui::IsMouseDown(ImGuiMouseButton_Left) && is_hovered && proximity_flag && !disable && !reshaping) {
				if (mouse_pos_in_canvas.x != old_mouse_pos.x || mouse_pos_in_canvas.y != old_mouse_pos.y) {
					points.push_back(mouse_pos_in_canvas);
					fragments.back().raw_data.push_back(mouse_pos_in_canvas);
				}
				drawing = true;
			}

			//actions taken on mouse release 
			if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
				//actions taken when finished drawing
				if (drawing == true) {
					points = LinearInterpolation(points, 1); //interpolate data such that it is evenly spaced
					fragments.back().raw_data = LinearInterpolation(fragments.back().raw_data, 1); //raw data is interpolated, but ONLY interpolated. 
					fragments.back().points = points; //update end fragment
					fragments.back().dir = Mixed; //needs to be processed to determine a direction
					new_state = true;
					drawing = false;
				}
				disable = false;
			}
			draw_list->PushClipRect(canvas_p0, canvas_p1, true);
			// Draw grid
			if (opt_enable_grid)
			{
				const float GRID_STEP = 64.0f;
				for (float x = fmodf(scrolling.x, GRID_STEP); x < canvas_sz.x; x += GRID_STEP)
					draw_list->AddLine(ImVec2(canvas_p0.x + x, canvas_p0.y), ImVec2(canvas_p0.x + x, canvas_p1.y), IM_COL32(200, 200, 200, 40));
				for (float y = fmodf(scrolling.y, GRID_STEP); y < canvas_sz.y; y += GRID_STEP)
					draw_list->AddLine(ImVec2(canvas_p0.x, canvas_p0.y + y), ImVec2(canvas_p1.x, canvas_p0.y + y), IM_COL32(200, 200, 200, 40));
			}
			//draw all the lines in the points vector
			//This draw section actually contains some data on how the drawn data should look, strictly based off of the user's inputs
			//Thus, the data is not stored in state, it is just modified here by how the user interacts with the program
			static ImU32 forward_color = IM_COL32(0, 255, 0, 255);
			static ImU32 forward_trans = IM_COL32(200, 100, 200, 150);
			static ImU32 backward_color = IM_COL32(255, 0, 0, 255);
			static ImU32 backward_trans = IM_COL32(100, 200, 200, 150);
			static ImU32 mixed_color = IM_COL32(255, 255, 0, 255);
			static ImU32 mixed_trans = IM_COL32(100, 100, 255, 150);
			float thickness = 2.0f;
			static float thck_modi = 3.5f;
			for (int i = 0; i < fragments.size(); i++) {
				if (i == clicked_frag)
					thickness = thickness * thck_modi;
				else
					thickness = 2.0f;
				if (fragments[i].dir & Forward && direction_colors) {
					for (int n = 0; (n + 1) < fragments[i].points.size(); n++)
						draw_list->AddLine(ImVec2(origin.x + fragments[i].points[n].x, origin.y + fragments[i].points[n].y), ImVec2(origin.x + fragments[i].points[n + 1].x, origin.y + fragments[i].points[n + 1].y), forward_color, thickness);
					//					for (int n = 0; (n + 1) < fragments[i].raw_data.size(); n+=2)
					//						draw_list->AddLine(ImVec2(origin.x + fragments[i].raw_data[n].x, origin.y + fragments[i].raw_data[n].y), ImVec2(origin.x + fragments[i].raw_data[n + 1].x, origin.y + fragments[i].raw_data[n + 1].y), forward_trans, thickness);
				}
				else if (fragments[i].dir & Backward && direction_colors) {
					for (int n = 0; (n + 1) < fragments[i].points.size(); n++)
						draw_list->AddLine(ImVec2(origin.x + fragments[i].points[n].x, origin.y + fragments[i].points[n].y), ImVec2(origin.x + fragments[i].points[n + 1].x, origin.y + fragments[i].points[n + 1].y), backward_color, thickness);
					//				for (int n = 0; (n + 1) < fragments[i].raw_data.size(); n+=2)
					//					draw_list->AddLine(ImVec2(origin.x + fragments[i].raw_data[n].x, origin.y + fragments[i].raw_data[n].y), ImVec2(origin.x + fragments[i].raw_data[n + 1].x, origin.y + fragments[i].raw_data[n + 1].y), backward_trans, thickness);
				}
				else {
					for (int n = 0; (n + 1) < fragments[i].points.size(); n++)
						draw_list->AddLine(ImVec2(origin.x + fragments[i].points[n].x, origin.y + fragments[i].points[n].y), ImVec2(origin.x + fragments[i].points[n + 1].x, origin.y + fragments[i].points[n + 1].y), mixed_color, thickness);
					//			for (int n = 0; (n + 1) < fragments[i].raw_data.size(); n+=2)
					//				draw_list->AddLine(ImVec2(origin.x + fragments[i].raw_data[n].x, origin.y + fragments[i].raw_data[n].y), ImVec2(origin.x + fragments[i].raw_data[n + 1].x, origin.y + fragments[i].raw_data[n + 1].y), mixed_trans, thickness);
				}
				vector<ImVec2> frag_points = fragments[i].points;
				if (!frag_points.empty())
					draw_list->AddCircle(ImVec2(origin.x + frag_points.back().x, origin.y + frag_points.back().y), 5, IM_COL32(0, 255, 255, 255), 20, 1.0);
			}



			draw_list->PopClipRect();
			draw_splitter.Merge(draw_list);
			ImGui::End();

		}
		//rendering
		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		glfwSwapBuffers(window);
	}

	//End Program
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

