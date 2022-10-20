#include <iostream>
#include <cmath>
#include <cstdio>
#include <random>
#include <fstream>
#include <algorithm>

// HDK/Imath stuff
#include <UT/UT_Vector2.h> // Ensure we are using HDK stuff
#include <Imath/ImathBoxAlgo.h> // This will include both ImathBox.h and ImathVec.h
#include <hboost/format.hpp>

using namespace std;

typedef std::vector<size_t> PointIndicesType;
typedef std::vector<Imath::V2f> PointArrayType;

PointArrayType globalPA;

static bool compareX(const size_t& p, const size_t& q) {
	return globalPA[p].x < globalPA[q].x;
}

static bool compareY(const size_t& p, const size_t& q) {
	return globalPA[p].y < globalPA[q].y;
}

typedef std::pair<size_t, size_t> SSPairType;
typedef std::pair<SSPairType, float> DResultType;
typedef std::vector<DResultType> DResultsContainer;

std::ostream & operator << (std::ostream &out, const DResultType &dresult)
{
	out << hboost::format("[%1%,%2%]:%3%") % dresult.first.first % dresult.first.second % dresult.second;
    return out;
}

std::ostream & operator << (std::ostream &out, const PointIndicesType &point_indices)
{
	out << "[ ";
	for (PointIndicesType::value_type p: point_indices)
		out << p << " ";
	out << "]";
    return out;
}

class MyClosestPair {
public:


	DResultType closestPair(PointIndicesType& L) {
		int n = L.size();
		if (n==0)
		{
			return DResultType
					{SSPairType{std::numeric_limits<SSPairType::first_type>::infinity(), std::numeric_limits<SSPairType::second_type>::infinity()}, std::numeric_limits<float>::infinity()};
		}
		printf("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ INITIAL L.size() = %ld\n",L.size());
		std::sort(L.begin(), L.end(), compareX);
		PointIndicesType res_points;
		DResultType d, d1, d2;
		printf("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ n = %d\n",n);
		int mid_point = n / 2;
		printf("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ mid_point = %d\n",mid_point);
		printf("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ L.size() = %ld\n",L.size());
		printf("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ L[mid_point] = %ld\n",L[mid_point]);
		printf("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ globalPA.size() = %ld\n",globalPA.size());
		float mid = globalPA[L[mid_point]][0];
		if (n <= 1) {
			return DResultType
					{SSPairType{L[0], L[0]}, std::numeric_limits<float>::infinity()};
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
			if (globalPA[L1[i]][0] > mid - d.second) {
				L_strip.push_back(L1[i]);
			}
		}
		for (int i = 0; i < L2.size(); i++) {
			if (globalPA[L2[i]][0] < mid + d.second) {
				L_strip.push_back(L2[i]);
			}
		}
		sort(L_strip.begin(), L_strip.end(), compareY);
		L = merge(L1, L2, false);
		if (L_strip.size() <= 1) {
			return d;
		} else {
			printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
			DResultType d_strip = DResultType
					{SSPairType{0, 0}, std::numeric_limits<float>::infinity()};
			for (int i = 0; i < L_strip.size(); i++) {
				// for (int j = i + 1; j < L_strip.size(); j++) {
				for (int j = i + 1; j < i + 8; j++) {
					if (j >= L_strip.size()) {
						break;
					}
					DResultType dist_successor = dist(L_strip[i], L_strip[j]);
					d_strip = dist_successor.second < d_strip.second ? dist_successor : d_strip;
				}
			}
			return d.second < d_strip.second ? d : d_strip;
		}
	}


	void closestPair(PointIndicesType& L,DResultsContainer& results) {
		DResultType ret = closestPair(L);
		std::cout << "ret : " << ret << std::endl;
		results.push_back(ret);
		if (std::isinf(ret.second))
		{
			std::cout << "ret : INFINITY" << std::endl;
		}
		if (L.size()<2)
			return;
		PointIndicesType L_san_ret;
		std::copy_if(L.begin(), L.end(),
		                 std::back_inserter(L_san_ret),
		                 [ret](size_t i) {
			bool status0 = (i != ret.first.first);
			bool status1 = (i != ret.first.second);
			bool status = status0 && status1;
			printf("status0 %s\n",(status0?"true":"false"));
			printf("status1 %s\n",(status1?"true":"false"));
			printf("status %s\n",(status?"true":"false"));
			return status;
		});
		printf("ret %ld %ld %f\n",ret.first.first,ret.first.second,ret.second);
		printf("closestPair : L -> %ld\n",L.size());
		printf("closestPair : L_san_ret -> %ld\n",L_san_ret.size());
		std::cout << "L : " << L << std::endl;
		std::cout << "L_san_ret : " << L_san_ret << std::endl;
		closestPair(L_san_ret,results);
	}

	void print_points(const PointIndicesType& points) {
		for (size_t p: points) {
			std::cout << globalPA[p];
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
			if (globalPA[L1[L1_ptr]][coor] <= globalPA[L2[L2_ptr]][coor]) {
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

	DResultType dist(const size_t& p, const size_t& q) {
		// float xDiff = p[0] - q[0];
		// float yDiff = p[1] - q[1];

		float xDiff = globalPA[p][0] - globalPA[q][0];
		float yDiff = globalPA[p][1] - globalPA[q][1];

		return DResultType
		{SSPairType {p, q}, sqrt(xDiff * xDiff + yDiff * yDiff)};
	}


}; // MyClosest Pair

void pointsFromFile(const std::string& filepath)
{
	PointIndicesType point_indices;
	std::cout << hboost::format("Processing using filepath = '%1%'") % filepath << std::endl;
	FILE *fp = fopen(filepath.c_str(),"r");
	if (fp)
	{
		char *response = NULL;
		size_t len;
		ssize_t status;
		size_t index = 0;
		while ((status = getline(&response, &len, fp)) > 1)
		{
			printf("status = '%ld'\n",status);
			printf("response = '%s'\n",response);
			float u,v;
			sscanf(response,"%f %f",&u,&v);
			// sscanf(response,"%f",&v);
			printf("UV = [%f,%f]\n",u,v);
			globalPA.push_back(Imath::V2f(u,v));
			point_indices.push_back(index++);
		}

		fclose(fp);


		MyClosestPair cp;
		DResultsContainer results;

		cp.closestPair(point_indices,results);
	}
}

void randomPointsToGnuPlot()
{
	int point_num = 100000;
	PointIndicesType point_indices(point_num);
	int rangeS = 0;
	int rangeE = 100000;
	srand(time(0));
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	std::uniform_real_distribution<> distr(rangeS, rangeE);
	// vector<Imath::V2f> points;
	clock_t tStart = clock();
	for (int i = 0; i != point_num; i++) {
		vector<float> coor{(float)distr(generator), (float)distr(generator)};
		// points.push_back(Imath::V2f(coor[0],coor[1]));
		globalPA.push_back(Imath::V2f(coor[0],coor[1]));
		point_indices[i] = i;
	}
	// printf("Time taken to initialize: \033[0;32m%fs\033[0m\n\n", (float)(clock() - tStart)/CLOCKS_PER_SEC);
	std::cout << hboost::format("Time taken to initialize: \033[0;32m%fs\033[0m\n\n") % ((float)(clock() - tStart)/CLOCKS_PER_SEC);
	// printf("\nThere are \033[1;31m%d\033[0m points, whose coordinates are in the interval \033[1;31m[%d, %d]\033[0m.\n\n", point_num, rangeS, rangeE);
	std::cout << hboost::format("\nThere are \033[1;31m%d\033[0m points, whose coordinates are in the interval \033[1;31m[%d, %d]\033[0m.\n\n") % point_num % rangeS % rangeE;
	tStart = clock();
	MyClosestPair cp;
	std::cout << hboost::format("\npi.size(BEFORE) = %1%") % point_indices.size() << std::endl;
	DResultType res = cp.closestPair(point_indices);
	std::cout << hboost::format("\npi.size(AFTER) = %1%") % point_indices.size() << std::endl;
	//printf("The closest pair of points is: \033[0;33m%s\033[0m and \033[0;33m%s\033[0m, with the distance \033[1;31m%f\033[0m.\n\n",
	//        res.first.first.toString().c_str(), res.first.second.toString().c_str(), res.second);
	std::cout << hboost::format("The closest pair of points is: \033[0;33m%s\033[0m and \033[0;33m%s\033[0m, with the distance \033[1;31m%f\033[0m.\n\n")
        		  % res.first.first % res.first.second % res.second;
	// printf("Time taken: \033[0;32m%fs\033[0m\n\n", (float)(clock() - tStart)/CLOCKS_PER_SEC);
	std::cout << hboost::format("Time taken: \033[0;32m%fs\033[0m\n\n") % ((float)(clock() - tStart)/CLOCKS_PER_SEC);
	std::ofstream file, fres;
	file.open("plot.txt");
	for (size_t p: point_indices) {
		file << globalPA[p][0] << "  " << globalPA[p][1] << std::endl;
	}
	file.close();
	fres.open("res.txt");
	fres << globalPA[res.first.first][0] << "  " << globalPA[res.first.first][1] << std::endl;
	fres << globalPA[res.first.second][0] << "  " << globalPA[res.first.second][1] << std::endl;
	fres.close();
	// HDK's Qt breaks gnuplot
	// system("gnuplot -e \"plot 'plot.txt' using 1:2 pt 7 ps 1 title 'Points', 'res.txt' using 1:2 pt 7 ps 1 lc rgb 'red' title 'Result'; pause -1\"");
}
/**
 * Run in terminal to compile:
 $ g++ -c merge_sort.cpp 2d_problem.cpp main.cpp && g++ merge_sort.o 2d_problem.o main.o -o my_program && ./my_program
 * */
int main(int argc, char** argv)
{
	if (argc==1)
		randomPointsToGnuPlot();
	else {
		pointsFromFile(argv[1]);
	}
	return 0;
}
