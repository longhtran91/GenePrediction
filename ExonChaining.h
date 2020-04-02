#pragma once
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include "Coordinate_Utility.h"
#include <QProgressDialog>

class ExonChaining
{
private:
	std::list<Interval_Coordinate> intervals;
    bool cancelled = false;

public:
	ExonChaining();
    ExonChaining(const std::list<Interval_Coordinate> &exons, QProgressDialog &p);
	std::list<Interval_Coordinate> get_intervals();
	~ExonChaining();
};
