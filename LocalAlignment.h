#pragma once
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <list>
#include <vector>
#include <unordered_set>
#include <array>
#include <stack>
#include "Coordinate_Utility.h"
#include <QProgressDialog>

#define DIAGONAL_PATH 2
#define VERTICAL_PATH 1
#define HORIZONTAL_PATH 4
#define MATCH_CHAR '|'
#define MISMATCH_CHAR '.'
#define SPACE_CHAR ' '
#define GAP_CHAR '-'

class LocalAlignment
{

private:
	std::list<Interval_Coordinate> exons;
	std::string dna_sequence;
	std::string dna_template;
	unsigned int  **trace_matrix;
    unsigned int row_size, col_size, num_row_executed;
    int match, mismatch, gap, threshold;

public:
	LocalAlignment();
    LocalAlignment(const std::string &dna_sequence, const std::string &dna_template, const int &match, const int &mismatch, const int &gap, const int &threshold, QProgressDialog &p);
    //LocalAlignment(const LocalAlignment &la);
	~LocalAlignment();
	std::list<Interval_Coordinate> getExons() const;
	std::list<std::array<std::string, 3>> print_alignment(const Coordinate &start, const Coordinate &end);	
};
