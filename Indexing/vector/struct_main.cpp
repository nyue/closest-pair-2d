#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <iterator>     // std::ostream_iterator

// HDK/Imath stuff
#include <UT/UT_Vector2.h> // Ensure we are using HDK stuff
#include <Imath/ImathBoxAlgo.h> // This will include both ImathBox.h and ImathVec.h
#include <hboost/format.hpp>
#include <hboost/bind.hpp>

using namespace std;

typedef std::vector<int> PointIndicesType;
typedef std::pair<int, int> IndicesPairType;

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
	int index() const { return _index;}
protected:
	int _index; // mimic OpenSubdiv refiner face indices which is a Vtr::Index of type `int`
	Imath::Box2f _bbox;
	Imath::V2f _centroid;
};
typedef std::vector<OSDFaceData> OSDFaceDataArrayType;

class MyClosestPair {
	const OSDFaceDataArrayType& _faceData;

	bool compareX(const size_t& p, const size_t& q) {
		return _faceData[p].centroid().x < _faceData[q].centroid().x;
	}

	bool compareY(const size_t& p, const size_t& q) {
		return _faceData[p].centroid().y < _faceData[q].centroid().y;
	}

public:

	typedef std::pair<IndicesPairType, float> DResult;

	MyClosestPair(const OSDFaceDataArrayType& data)
	:_faceData(data){

	}
	DResult closestPair(PointIndicesType& L) {
		std::sort(L.begin(), L.end(), hboost::bind(&MyClosestPair::compareX, this, _1, _2));
		PointIndicesType res_points;
		DResult d, d1, d2;
		int n = L.size();
		int mid_point = n / 2;
		float mid = _faceData[L[mid_point]].centroid()[0];
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
			if (_faceData[L1[i]].centroid()[0] > mid - d.second) {
				L_strip.push_back(L1[i]);
			}
		}
		for (int i = 0; i < L2.size(); i++) {
			if (_faceData[L2[i]].centroid()[0] < mid + d.second) {
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
			std::cout << _faceData[p].centroid();
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
			if (_faceData[L1[L1_ptr]].centroid()[coor] <= _faceData[L2[L2_ptr]].centroid()[coor]) {
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

		float xDiff = _faceData[p].centroid()[0] - _faceData[q].centroid()[0];
		float yDiff = _faceData[p].centroid()[1] - _faceData[q].centroid()[1];

		return std::pair<std::pair<size_t, size_t>, float>
		{std::pair<size_t, size_t> {p, q}, sqrt(xDiff * xDiff + yDiff * yDiff)};
	}


}; // MyClosest Pair

void generate_face_data(OSDFaceDataArrayType& faceData)
{
	int point_num = 10000;
	int rangeS = 0;
	int rangeE = 1;
	srand(time(0));
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	std::uniform_real_distribution<> distr(rangeS, rangeE);
	clock_t tStart = clock();
	for (int i = 0; i != point_num; i++) {
		faceData.push_back(OSDFaceData());
		faceData.back().init(i, Imath::Box2f(Imath::V2f((float)distr(generator), (float)distr(generator)),Imath::V2f((float)distr(generator), (float)distr(generator))));
	}
	std::cout << hboost::format("Time taken to initialize: \033[0;32m%fs\033[0m\n\n") % ((float)(clock() - tStart)/CLOCKS_PER_SEC);
	std::cout << hboost::format("\nThere are \033[1;31m%d\033[0m bboxes, whose coordinates are in the interval \033[1;31m[%d, %d]\033[0m.\n\n") % point_num % rangeS % rangeE;
}

void get_point_indices_from_face_data(const OSDFaceDataArrayType& faceData,
		PointIndicesType& o_pi)
{
	for (auto item:faceData)
	{
		o_pi.push_back(item.index());
	}
}

void remove_processed_pair(const IndicesPairType& processedPair, PointIndicesType& indicesData)
{
	indicesData.erase(std::remove(indicesData.begin(), indicesData.end(), processedPair.first), indicesData.end());
	indicesData.erase(std::remove(indicesData.begin(), indicesData.end(), processedPair.second), indicesData.end());
}

/**
 * Run in terminal to compile:
 $ g++ -c merge_sort.cpp 2d_problem.cpp main.cpp && g++ merge_sort.o 2d_problem.o main.o -o my_program && ./my_program
 * */
int main() {

#if TEST_REMOVE_PAIR
	PointIndicesType indicesData;
	indicesData.push_back(3);
	indicesData.push_back(1);
	indicesData.push_back(4);
	indicesData.push_back(5);
	indicesData.push_back(9);
	indicesData.push_back(2);
	indicesData.push_back(6);

	std::copy(indicesData.begin(), indicesData.end(), std::ostream_iterator<PointIndicesType::value_type>(std::cout, " "));
	std::cout << std::endl;

	IndicesPairType processed_pair(1,9);

	remove_processed_pair(processed_pair, indicesData);

	std::copy(indicesData.begin(), indicesData.end(), std::ostream_iterator<PointIndicesType::value_type>(std::cout, " "));
	std::cout << std::endl;
#else // TEST_REMOVE_PAIR
	OSDFaceDataArrayType faceData;
	PointIndicesType indicesData;

	generate_face_data(faceData);
	get_point_indices_from_face_data(faceData,indicesData);
	clock_t tStart = clock();
	MyClosestPair cp(faceData);
	while (indicesData.size()>1) {
		// std::cout << hboost::format("\npi.size(BEFORE) = %1%") % indicesData.size() << std::endl;
		MyClosestPair::DResult res = cp.closestPair(indicesData);
		remove_processed_pair(res.first,indicesData);
		// std::cout << hboost::format("\npi.size(AFTER) = %1%") % indicesData.size() << std::endl;
		// std::cout << hboost::format("The closest pair of points is: \033[0;33m%s\033[0m and \033[0;33m%s\033[0m, with the distance \033[1;31m%f\033[0m.\n\n")
	    //    		  % res.first.first % res.first.second % res.second;
	}
	std::cout << hboost::format("Time taken: \033[0;32m%fs\033[0m\n\n") % ((float)(clock() - tStart)/CLOCKS_PER_SEC);
#ifdef ENABLE_GNUPLOT
	std::ofstream file, fres;
	file.open("plot_struct.txt");
	for (size_t p: indicesData) {
		file << faceData[p].centroid()[0] << "  " << faceData[p].centroid()[1] << std::endl;
	}
	file.close();
	fres.open("res_struct.txt");
	fres << faceData[res.first.first].centroid()[0] << "  " << faceData[res.first.first].centroid()[1] << std::endl;
	fres << faceData[res.first.second].centroid()[0] << "  " << faceData[res.first.second].centroid()[1] << std::endl;
	fres.close();
	// HDK's Qt breaks gnuplot
	// system("gnuplot -e \"plot 'plot.txt' using 1:2 pt 7 ps 1 title 'Points', 'res.txt' using 1:2 pt 7 ps 1 lc rgb 'red' title 'Result'; pause -1\"");
#endif // ENABLE_GNUPLOT
#endif // TEST_REMOVE_PAIR

	return 0;
}
