#include <iostream>
#include <cmath>
#include <random>
// #include "../include/2d_problem.h"
#include <fstream>
// #include "matplotlibcpp.h"
#include <algorithm>

using namespace std;

// point stuff START
class Point {
  public:
    // only index 0 and 1 of the vector for the 2d problem
    std::vector<double> coordinate;
    Point(std::vector<double> coordinate) {
      this->coordinate = coordinate;
    }
    // Point(double x, double y) {
    //   coordinate.push_back(x);
    //   coordinate.push_back(y);
    // }
    Point(){}
    std::string toString() {
      return "(" +  std::to_string(coordinate[0]) + "," + std::to_string(coordinate[1]) + ")";
    }
    friend bool operator== (const Point& p1, const Point& p2);
    friend std::ostream& operator<<(std::ostream &strm, const Point &a) {
      return strm << "(" << a.coordinate[0] << "," << a.coordinate[1] << ")";
    }
};

inline bool operator== (const Point& p1, const Point& p2) {
  return p1.coordinate == p2.coordinate;
}

// point stuff END

double compareX(Point p, Point q) {
	return p.coordinate[0] < q.coordinate[0];
}

double compareY(Point p, Point q) {
	return p.coordinate[1] < q.coordinate[1];
}


class MyClosestPair {

public:
	// merge_sort stuff START
	/*
	vector<Point> merge_sort(vector<Point>& L, bool x);
	vector<Point> merge(vector<Point> const &L1, vector<Point> const &L2, bool x);
	 */
	vector<Point> merge_sort(vector<Point> &L, bool x) {
		if (L.size() <= 1) {
			return L;
		}

		vector<Point> L1, L2;
		int mid = ceil(L.size() / 2);

		for (auto i = L.begin(); i != L.begin() + mid; i++) {
			L1.push_back(*i);
		}

		for (auto i = L.begin() + mid; i != L.end(); i++) {
			L2.push_back(*i);
		}

		L1 = merge_sort(L1, x);
		L2 = merge_sort(L2, x);
		L = merge(L1, L2, x);
		return L;
	}

	vector<Point> merge(vector<Point> const &L1, vector<Point> const &L2, bool x) {
		int coor = x ? 0 : 1;
		vector<Point> res;
		int L1_ptr = 0;
		int L2_ptr = 0;
		while (L1_ptr < L1.size() and L2_ptr < L2.size()) {
			if (L1[L1_ptr].coordinate[coor] <= L2[L2_ptr].coordinate[coor]) {
				res.push_back(L1[L1_ptr++]);
			} else {
				res.push_back(L2[L2_ptr++]);
			}
		}
		if (L1_ptr == L1.size()) {
			res.insert(res.end(), L2.begin() + L2_ptr, L2.end());
		} else if (L2_ptr == L2.size()) {
			res.insert(res.end(), L1.begin() + L1_ptr, L1.end());
		}
		return res;
	}

	// merge_sort stuff END

	// 2d_problem stuff START
	typedef std::pair<std::pair<Point, Point>, double> DResult;

	/*
		DResult closestPair(vector<Point>& L);
		void print_points(vector<Point> points);
	 */

	DResult dist(Point p, Point q) {
		double xDiff = p.coordinate[0] - q.coordinate[0];
		double yDiff = p.coordinate[1] - q.coordinate[1];
		return std::pair<std::pair<Point, Point>, double>
		{std::pair<Point, Point> {p, q}, sqrt(xDiff * xDiff + yDiff * yDiff)};
	}
	/*
	void print_points(vector<Point> points);
	*/

	DResult closestPair(vector<Point>& L) {
		std::sort(L.begin(), L.end(), compareX);
		vector<Point> res_points;
		DResult d, d1, d2;
		int n = L.size();
		int mid_point = n / 2;
		double mid = L[mid_point].coordinate[0];
		if (n <= 1) {
			return DResult
					{std::pair<Point, Point>{Point(), Point()}, std::numeric_limits<double>::infinity()};
		} else if (n == 2) {
			return dist(L[0], L[1]);
		}
		// merge_sort(L, true);
		vector<Point> L1, L2, L_strip;
		for (int i = 0; i < n; i++) {
			if (i < mid_point) {
				L1.push_back(L[i]);
			} else {
				L2.push_back(L[i]);
			}
		}
		d1 = closestPair(L1);
		d2 = closestPair(L2);
		d = d1.second < d2.second ? d1 : d2;

		for (int i = L1.size() - 1; i >= 0; i--) {
			if (L1[i].coordinate[0] > mid - d.second) {
				L_strip.push_back(L1[i]);
			}
		}
		for (int i = 0; i < L2.size(); i++) {
			if (L2[i].coordinate[0] < mid + d.second) {
				L_strip.push_back(L2[i]);
			}
		}
		sort(L_strip.begin(), L_strip.end(), compareY);
		L = merge(L1, L2, false);
		if (L_strip.size() <= 1) {
			return d;
		} else {
			DResult d_strip = DResult
					{std::pair<Point, Point>{Point(), Point()}, std::numeric_limits<double>::infinity()};
			for (int i = 0; i < L_strip.size(); i++) {
				// for (int j = i + 1; j < L_strip.size(); j++) {
				for (int j = i + 1; j < i + 8; j++) {
					if (j >= L_strip.size()) {
						break;
					}
					DResult dist_successor = dist(L_strip[i], L_strip[j]);
					d_strip = dist_successor.second < d_strip.second ? dist_successor : d_strip;
				}
			}
			return d.second < d_strip.second ? d : d_strip;
		}
	}


	void print_points(vector<Point> points) {
		for (Point p: points) {
			std::cout << p;
		}
		std::cout << endl;
	}
}; // MyClosest Pair

// 2d_problem stuff END

/**
 * Run in terminal to compile:
 $ g++ -c merge_sort.cpp 2d_problem.cpp main.cpp && g++ merge_sort.o 2d_problem.o main.o -o my_program && ./my_program 
 * */
int main() {
  int point_num = 100000;
  int rangeS = 0;
  int rangeE = 100000;
  srand(time(0));
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<> distr(rangeS, rangeE);
  vector<Point> points;
  clock_t tStart = clock();
  for (int i = 0; i != point_num; i++) {
    vector<double> coor{distr(generator), distr(generator)};
    points.push_back(Point(coor));
  }
  printf("Time taken to initialize: \033[0;32m%fs\033[0m\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  printf("\nThere are \033[1;31m%d\033[0m points, whose coordinates are in the interval \033[1;31m[%d, %d]\033[0m.\n\n", point_num, rangeS, rangeE);
  tStart = clock();
  MyClosestPair cp;
  MyClosestPair::DResult res = cp.closestPair(points);
  printf("The closest pair of points is: \033[0;33m%s\033[0m and \033[0;33m%s\033[0m, with the distance \033[1;31m%f\033[0m.\n\n",
          res.first.first.toString().c_str(), res.first.second.toString().c_str(), res.second);
  printf("Time taken: \033[0;32m%fs\033[0m\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  std::ofstream file, fres;
  file.open("plot.txt");
  for (Point p: points) {
    file << p.coordinate[0] << "  " << p.coordinate[1] << std::endl;
  }
  file.close();
  fres.open("res.txt"); 
  fres << res.first.first.coordinate[0] << "  " << res.first.first.coordinate[1] << std::endl;
  fres << res.first.second.coordinate[0] << "  " << res.first.second.coordinate[1] << std::endl;
  fres.close();
  system("gnuplot -e \"plot 'plot.txt' using 1:2 pt 7 ps 1 title 'Points', 'res.txt' using 1:2 pt 7 ps 1 lc rgb 'red' title 'Result'; pause -1\"");
  return 0;
}
