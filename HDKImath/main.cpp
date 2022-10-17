#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>

// HDK/Imath stuff
#include <UT/UT_Vector2.h> // Ensure we are using HDK stuff
#include <Imath/ImathBoxAlgo.h> // This will include both ImathBox.h and ImathVec.h
#include <hboost/format.hpp>

using namespace std;

bool compareX(const Imath::V2f& p, const Imath::V2f& q) {
	return p[0] < q[0];
}

bool compareY(const Imath::V2f& p, const Imath::V2f& q) {
	return p[1] < q[1];
}


class MyClosestPair {
public:

	typedef std::pair<std::pair<Imath::V2f, Imath::V2f>, float> DResult;

	DResult closestPair(vector<Imath::V2f>& L) {
		std::sort(L.begin(), L.end(), compareX);
		vector<Imath::V2f> res_points;
		DResult d, d1, d2;
		int n = L.size();
		int mid_point = n / 2;
		float mid = L[mid_point][0];
		if (n <= 1) {
			return DResult
					{std::pair<Imath::V2f, Imath::V2f>{Imath::V2f(), Imath::V2f()}, std::numeric_limits<float>::infinity()};
		} else if (n == 2) {
			return dist(L[0], L[1]);
		}
		// merge_sort(L, true);
		vector<Imath::V2f> L1, L2, L_strip;
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
			if (L1[i][0] > mid - d.second) {
				L_strip.push_back(L1[i]);
			}
		}
		for (int i = 0; i < L2.size(); i++) {
			if (L2[i][0] < mid + d.second) {
				L_strip.push_back(L2[i]);
			}
		}
		sort(L_strip.begin(), L_strip.end(), compareY);
		L = merge(L1, L2, false);
		if (L_strip.size() <= 1) {
			return d;
		} else {
			DResult d_strip = DResult
					{std::pair<Imath::V2f, Imath::V2f>{Imath::V2f(), Imath::V2f()}, std::numeric_limits<float>::infinity()};
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


	void print_points(const vector<Imath::V2f>& points) {
		for (Imath::V2f p: points) {
			std::cout << p;
		}
		std::cout << endl;
	}
private:

	vector<Imath::V2f> merge_sort(vector<Imath::V2f> &L, bool x) {
		if (L.size() <= 1) {
			return L;
		}

		vector<Imath::V2f> L1, L2;
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

	vector<Imath::V2f> merge(vector<Imath::V2f> const &L1, vector<Imath::V2f> const &L2, bool x) {
		int coor = x ? 0 : 1;
		vector<Imath::V2f> res;
		int L1_ptr = 0;
		int L2_ptr = 0;
		while (L1_ptr < L1.size() and L2_ptr < L2.size()) {
			if (L1[L1_ptr][coor] <= L2[L2_ptr][coor]) {
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

	DResult dist(const Imath::V2f& p, const Imath::V2f& q) {
		float xDiff = p[0] - q[0];
		float yDiff = p[1] - q[1];
		return std::pair<std::pair<Imath::V2f, Imath::V2f>, float>
		{std::pair<Imath::V2f, Imath::V2f> {p, q}, sqrt(xDiff * xDiff + yDiff * yDiff)};
	}


}; // MyClosest Pair

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
  vector<Imath::V2f> points;
  clock_t tStart = clock();
  for (int i = 0; i != point_num; i++) {
    vector<float> coor{(float)distr(generator), (float)distr(generator)};
    points.push_back(Imath::V2f(coor[0],coor[1]));
  }
  // printf("Time taken to initialize: \033[0;32m%fs\033[0m\n\n", (float)(clock() - tStart)/CLOCKS_PER_SEC);
  std::cout << hboost::format("Time taken to initialize: \033[0;32m%fs\033[0m\n\n") % ((float)(clock() - tStart)/CLOCKS_PER_SEC);
  // printf("\nThere are \033[1;31m%d\033[0m points, whose coordinates are in the interval \033[1;31m[%d, %d]\033[0m.\n\n", point_num, rangeS, rangeE);
  std::cout << hboost::format("\nThere are \033[1;31m%d\033[0m points, whose coordinates are in the interval \033[1;31m[%d, %d]\033[0m.\n\n") % point_num % rangeS % rangeE;
  tStart = clock();
  MyClosestPair cp;
  MyClosestPair::DResult res = cp.closestPair(points);
  //printf("The closest pair of points is: \033[0;33m%s\033[0m and \033[0;33m%s\033[0m, with the distance \033[1;31m%f\033[0m.\n\n",
  //        res.first.first.toString().c_str(), res.first.second.toString().c_str(), res.second);
  std::cout << hboost::format("The closest pair of points is: \033[0;33m%s\033[0m and \033[0;33m%s\033[0m, with the distance \033[1;31m%f\033[0m.\n\n")
          % res.first.first % res.first.second % res.second;
  // printf("Time taken: \033[0;32m%fs\033[0m\n\n", (float)(clock() - tStart)/CLOCKS_PER_SEC);
  std::cout << hboost::format("Time taken: \033[0;32m%fs\033[0m\n\n") % ((float)(clock() - tStart)/CLOCKS_PER_SEC);
  std::ofstream file, fres;
  file.open("plot.txt");
  for (Imath::V2f p: points) {
    file << p[0] << "  " << p[1] << std::endl;
  }
  file.close();
  fres.open("res.txt"); 
  fres << res.first.first[0] << "  " << res.first.first[1] << std::endl;
  fres << res.first.second[0] << "  " << res.first.second[1] << std::endl;
  fres.close();
  // HDK's Qt breaks gnuplot
  // system("gnuplot -e \"plot 'plot.txt' using 1:2 pt 7 ps 1 title 'Points', 'res.txt' using 1:2 pt 7 ps 1 lc rgb 'red' title 'Result'; pause -1\"");
  return 0;
}
