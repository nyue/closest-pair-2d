#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>

// HDK/Imath stuff
#include <UT/UT_Vector2.h> // Ensure we are using HDK stuff
#include <Imath/ImathBoxAlgo.h> // This will include both ImathBox.h and ImathVec.h
#include <hboost/format.hpp>
#include <hboost/bind.hpp>

using namespace std;

typedef std::vector<size_t> PointIndicesType;
class OSDFaceData {
public:
	OSDFaceData()
	:_index(-1)
	{
	}
	void init(int face_index, const Imath::Box2f& face_bbox)
	{
		_index = face_index;
		_bbox = face_bbox;
		_centroid = 0.5f*(_bbox.max-_bbox.min);
	}
	const Imath::V2f& centroid() const {
		assert(_index>=0);
		return _centroid;
	}
protected:
	int _index; // mimic OpenSubdiv refiner face indices which is a Vtr::Index of type `int`
	Imath::Box2f _bbox;
	Imath::V2f _centroid;
};
typedef std::vector<OSDFaceData> OSDFaceDataArrayType;

OSDFaceDataArrayType globalPA;

class MyClosestPair {
	bool compareX(const size_t& p, const size_t& q) {
		return globalPA[p].centroid().x < globalPA[q].centroid().x;
	}

	bool compareY(const size_t& p, const size_t& q) {
		return globalPA[p].centroid().y < globalPA[q].centroid().y;
	}

public:

	typedef std::pair<std::pair<size_t, size_t>, float> DResult;

	DResult closestPair(PointIndicesType& L) {
		std::sort(L.begin(), L.end(), hboost::bind(&MyClosestPair::compareX, this, _1, _2));
		PointIndicesType res_points;
		DResult d, d1, d2;
		int n = L.size();
		int mid_point = n / 2;
		float mid = globalPA[L[mid_point]].centroid()[0];
		if (n <= 1) {
			return DResult
					{std::pair<size_t, size_t>{0, 0}, std::numeric_limits<float>::infinity()};
		} else if (n == 2) {
			return dist(L[0], L[1]);
		}

		PointIndicesType L1, L2, L_strip;
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
			if (globalPA[L1[i]].centroid()[0] > mid - d.second) {
				L_strip.push_back(L1[i]);
			}
		}
		for (int i = 0; i < L2.size(); i++) {
			if (globalPA[L2[i]].centroid()[0] < mid + d.second) {
				L_strip.push_back(L2[i]);
			}
		}
		sort(L_strip.begin(), L_strip.end(), hboost::bind(&MyClosestPair::compareY, this, _1, _2));
		L = merge(L1, L2, false);
		if (L_strip.size() <= 1) {
			return d;
		} else {
			DResult d_strip = DResult
					{std::pair<size_t, size_t>{0, 0}, std::numeric_limits<float>::infinity()};
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


	void print_points(const PointIndicesType& points) {
		for (size_t p: points) {
			std::cout << globalPA[p].centroid();
		}
		std::cout << endl;
	}

private:

	PointIndicesType merge_sort(PointIndicesType&L, bool x) {
		if (L.size() <= 1) {
			return L;
		}

		PointIndicesType L1, L2;
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

	PointIndicesType merge(PointIndicesType const &L1, PointIndicesType const &L2, bool x) {
		int coor = x ? 0 : 1;
		PointIndicesType res;
		int L1_ptr = 0;
		int L2_ptr = 0;
		while (L1_ptr < L1.size() and L2_ptr < L2.size()) {
			if (globalPA[L1[L1_ptr]].centroid()[coor] <= globalPA[L2[L2_ptr]].centroid()[coor]) {
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

	DResult dist(const size_t& p, const size_t& q) {
		// float xDiff = p[0] - q[0];
		// float yDiff = p[1] - q[1];

		float xDiff = globalPA[p].centroid()[0] - globalPA[q].centroid()[0];
		float yDiff = globalPA[p].centroid()[1] - globalPA[q].centroid()[1];

		return std::pair<std::pair<size_t, size_t>, float>
		{std::pair<size_t, size_t> {p, q}, sqrt(xDiff * xDiff + yDiff * yDiff)};
	}


}; // MyClosest Pair

/**
 * Run in terminal to compile:
 $ g++ -c merge_sort.cpp 2d_problem.cpp main.cpp && g++ merge_sort.o 2d_problem.o main.o -o my_program && ./my_program 
 * */
int main() {
	int point_num = 100000;
	PointIndicesType pi(point_num);
	int rangeS = 0;
	int rangeE = 1;
	srand(time(0));
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	std::uniform_real_distribution<> distr(rangeS, rangeE);
	// vector<Imath::V2f> points;
	clock_t tStart = clock();
	for (int i = 0; i != point_num; i++) {
		// vector<float> coor{(float)distr(generator), (float)distr(generator)};
		// points.push_back(Imath::V2f(coor[0],coor[1]));
		// globalPA.push_back(Imath::V2f(coor[0],coor[1]));
		globalPA.push_back(OSDFaceData());
		globalPA.back().init(i, Imath::Box2f(Imath::V2f((float)distr(generator), (float)distr(generator)),Imath::V2f((float)distr(generator), (float)distr(generator))));
		pi[i] = i;
	}
	// printf("Time taken to initialize: \033[0;32m%fs\033[0m\n\n", (float)(clock() - tStart)/CLOCKS_PER_SEC);
	std::cout << hboost::format("Time taken to initialize: \033[0;32m%fs\033[0m\n\n") % ((float)(clock() - tStart)/CLOCKS_PER_SEC);
	// printf("\nThere are \033[1;31m%d\033[0m points, whose coordinates are in the interval \033[1;31m[%d, %d]\033[0m.\n\n", point_num, rangeS, rangeE);
	std::cout << hboost::format("\nThere are \033[1;31m%d\033[0m points, whose coordinates are in the interval \033[1;31m[%d, %d]\033[0m.\n\n") % point_num % rangeS % rangeE;
	tStart = clock();
	MyClosestPair cp;
	std::cout << hboost::format("\npi.size(BEFORE) = %1%") % pi.size() << std::endl;
	MyClosestPair::DResult res = cp.closestPair(pi);
	std::cout << hboost::format("\npi.size(AFTER) = %1%") % pi.size() << std::endl;
	//printf("The closest pair of points is: \033[0;33m%s\033[0m and \033[0;33m%s\033[0m, with the distance \033[1;31m%f\033[0m.\n\n",
	//        res.first.first.toString().c_str(), res.first.second.toString().c_str(), res.second);
	std::cout << hboost::format("The closest pair of points is: \033[0;33m%s\033[0m and \033[0;33m%s\033[0m, with the distance \033[1;31m%f\033[0m.\n\n")
        		  % res.first.first % res.first.second % res.second;
	// printf("Time taken: \033[0;32m%fs\033[0m\n\n", (float)(clock() - tStart)/CLOCKS_PER_SEC);
	std::cout << hboost::format("Time taken: \033[0;32m%fs\033[0m\n\n") % ((float)(clock() - tStart)/CLOCKS_PER_SEC);
	std::ofstream file, fres;
	file.open("plot_struct.txt");
	for (size_t p: pi) {
		file << globalPA[p].centroid()[0] << "  " << globalPA[p].centroid()[1] << std::endl;
	}
	file.close();
	fres.open("res_struct.txt");
	fres << globalPA[res.first.first].centroid()[0] << "  " << globalPA[res.first.first].centroid()[1] << std::endl;
	fres << globalPA[res.first.second].centroid()[0] << "  " << globalPA[res.first.second].centroid()[1] << std::endl;
	fres.close();
	// HDK's Qt breaks gnuplot
	// system("gnuplot -e \"plot 'plot.txt' using 1:2 pt 7 ps 1 title 'Points', 'res.txt' using 1:2 pt 7 ps 1 lc rgb 'red' title 'Result'; pause -1\"");
	return 0;
}
