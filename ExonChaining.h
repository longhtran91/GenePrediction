#pragma once
#include <array>
//#include <algorithm>
#include <map>
#include <set>
#include <unordered_set>
#include <stack>


using namespace std;

class ExonChaining
{
private:
	stack<array<unsigned int, 3>> intervals;

public:
	ExonChaining();
	ExonChaining(const unordered_set<array<unsigned int, 3>> &exons);
	stack<array<unsigned int, 3>> get_intervals();
	~ExonChaining();
	template <typename T>
	T find_max_itr(T start, T end);
};
