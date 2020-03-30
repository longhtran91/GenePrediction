#pragma once
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include "Coordinate_Utility.h"

class ExonChaining
{
private:
	std::list<Interval_Coordinate> intervals;

public:
	ExonChaining();
	ExonChaining(const std::list<Interval_Coordinate> &exons);
	std::list<Interval_Coordinate> get_intervals();
	~ExonChaining();
};
