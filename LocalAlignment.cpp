#include "LocalAlignment.h"
#include <chrono>
using namespace std;

LocalAlignment::LocalAlignment()
{

}

LocalAlignment::LocalAlignment(const string &dna_template, const string &dna_sequence, const int &match, const int &mismatch, const int &gap, const int &threshold, QProgressDialog &p)
{
    //this->read_config("var.config");
	this->dna_template = dna_template;
	this->dna_sequence = dna_sequence;
    this->match = match;
    this->mismatch = mismatch;
    this->gap = gap;
    this->threshold = threshold;
	this->row_size = dna_template.length() + 1;
	this->col_size = dna_sequence.length() + 1;
	unsigned int **scores_matrix = new unsigned int *[this->row_size];
    unordered_set<Coordinate> **origins_matrix = new unordered_set<Coordinate> *[this->row_size];
	this->trace_matrix = new unsigned int *[this->row_size];   


    unsigned int numTask = 0;
    bool cancelled = false;
	//fill the rest of matrix by dynamic programming local alignment
	for (unsigned int row = 0; row < this->row_size; ++row)
	{
        if (cancelled)
        {
            delete[] scores_matrix;
            delete[] origins_matrix;
            break;
        }
		//auto start = chrono::high_resolution_clock::now();
		scores_matrix[row] = new unsigned int[this->col_size];
        origins_matrix[row] = new unordered_set<Coordinate>[this->col_size];
		trace_matrix[row] = new unsigned int[this->col_size];
        num_row_executed = row;

		//rolling delete used/obsolete rows
		if (row >= 2)
		{
			delete[] scores_matrix[row - 2];
			delete[] origins_matrix[row - 2];
		}

		for (unsigned int col = 0; col < this->col_size; ++col)
		{
            if (p.wasCanceled())
            {
                cancelled = true;
                delete[] scores_matrix[row];
                delete[] scores_matrix[row-1];
                delete[] origins_matrix[row];
                delete[] origins_matrix[row-1];
                break;
            }

			//initialize first row and column to be 0
			if (row == 0 || col == 0)
			{
				scores_matrix[row][col] = 0;
				trace_matrix[row][col] = NULL;
				continue;
			}

			int vertical_score, horizontal_score, diagonal_score;

			horizontal_score = scores_matrix[row][col - 1] + this->gap;
			vertical_score = scores_matrix[row - 1][col] + this->gap;

			if (dna_sequence[col - 1] == dna_template[row - 1])	diagonal_score = scores_matrix[row - 1][col - 1] + this->match;
			else diagonal_score = scores_matrix[row - 1][col - 1] + this->mismatch;

			unsigned int score = max({ 0, vertical_score, horizontal_score, diagonal_score });
			scores_matrix[row][col] = score;

			trace_matrix[row][col] = 0;
			if (score == vertical_score)
			{
				trace_matrix[row][col] += VERTICAL_PATH;
                if (origins_matrix[row - 1][col].empty()) origins_matrix[row][col].insert(Coordinate(row, col));
                else origins_matrix[row][col].insert(origins_matrix[row - 1][col].cbegin(), origins_matrix[row - 1][col].cend());
			}
			if (score == diagonal_score)
			{
				trace_matrix[row][col] += DIAGONAL_PATH;
                if (origins_matrix[row - 1][col - 1].empty()) origins_matrix[row][col].insert(Coordinate(row, col));
                else origins_matrix[row][col].insert(origins_matrix[row - 1][col - 1].cbegin(), origins_matrix[row - 1][col - 1].cend());
			}
			if (score == horizontal_score)
			{
				trace_matrix[row][col] += HORIZONTAL_PATH;
                if (origins_matrix[row][col - 1].empty()) origins_matrix[row][col].insert(Coordinate(row, col));
                else origins_matrix[row][col].insert(origins_matrix[row][col - 1].cbegin(), origins_matrix[row][col - 1].cend());
            }

			if (score >= this->threshold)
			{
				for (const auto& elem : origins_matrix[row][col])
				{
					this->exons.push_back(Interval_Coordinate(Coordinate(elem.row, elem.col), Coordinate(row, col), score));
				}
			}
            p.setValue(numTask);
            ++numTask;
		}		
	}

    if (!cancelled)
    {
        for (int k = 1; k < 3; ++k)
        {
            delete[] scores_matrix[this->row_size - k];
            delete[] origins_matrix[this->row_size - k];
        }
        delete[] scores_matrix;
        delete[] origins_matrix;
    }

}
LocalAlignment::~LocalAlignment()
{
    for (unsigned int row = 0; row <= num_row_executed; ++row)
	{
		delete[] this->trace_matrix[row];
	}
	delete[] this->trace_matrix;
}

list<Interval_Coordinate> LocalAlignment::getExons() const
{
	return this->exons;
}

list<array<string, 3>> LocalAlignment::print_alignment(const Coordinate &start, const Coordinate &end)
{
	list<array<string, 3>> result;
	string ali_template;
	string ali_sequence;
	string mid_line;
	int index = 0;
	unsigned int trace = this->trace_matrix[end.row][end.col];
	if ((trace - trace % 2) % 4 == 2) this->trace_matrix[end.row][end.col] = 2;
	if (trace == NULL)	return result;

	if (end == start)
	{
		{
			char d = (this->dna_template[end.row - 1] == this->dna_sequence[end.col - 1] ? MATCH_CHAR : MISMATCH_CHAR);
			result.push_back({ string(1, this->dna_template[end.row - 1]), string(1, d),string(1, this->dna_sequence[end.col - 1]) });
			return result;
		}
	}

	stack<Coordinate> s;
	stack<int> si;
	si.push(index);
	s.push(end);
	Coordinate curr;
	Coordinate pre_curr;
	Coordinate tmp_curr;
	while (!s.empty())
	{
		bool col_gap = false, dead_end = false;
		curr = s.top();
		pre_curr = s.top();
		index = si.top();
		trace = this->trace_matrix[curr.row][curr.col];
		si.pop();
		s.pop();

		while (trace != NULL && pre_curr != start && !dead_end)
		{
			tmp_curr = curr;
			int size = ali_template.size();
			ali_template = ali_template.substr(size - index, size);
			ali_sequence = ali_sequence.substr(size - index, size);
            mid_line = mid_line.substr(size - index, size);

			if (trace % 2 == VERTICAL_PATH)
			{
				trace -= 1;
				col_gap = true;
				ali_template.insert(0, string(1, this->dna_template[curr.row - 1]));
				ali_sequence.insert(0, string(1, GAP_CHAR));
				mid_line.insert(0, string(1, SPACE_CHAR));
				--curr.row;
			}
			else if (trace % 4 == DIAGONAL_PATH)
			{
				trace -= 2;
				col_gap = false;
				char d = (this->dna_template[curr.row - 1] == this->dna_sequence[curr.col - 1] ? MATCH_CHAR : MISMATCH_CHAR);
				ali_template.insert(0, string(1, this->dna_template[curr.row - 1]));
				ali_sequence.insert(0, string(1, this->dna_sequence[curr.col - 1]));
				mid_line.insert(0, string(1, d));
				--curr.row;
				--curr.col;
			}
			else if (trace == HORIZONTAL_PATH)
			{
				if (col_gap) dead_end = true;
				else
				{
					trace -= 4;
					col_gap = false;
					ali_template.insert(0, string(1, GAP_CHAR));
					ali_sequence.insert(0, string(1, this->dna_sequence[curr.col - 1]));
					mid_line.insert(0, string(1, SPACE_CHAR));
					--curr.col;
				}
			}

			if (trace > 0 && !dead_end)
			{
				this->trace_matrix[tmp_curr.row][tmp_curr.col] = trace;
				s.push(tmp_curr);
				si.push(index);
			}
			pre_curr = tmp_curr;
			trace = this->trace_matrix[curr.row][curr.col];
			++index;
		}

		if (pre_curr == start)
		{
			result.push_back({ ali_template, mid_line, ali_sequence });
		}
	}
	return result;
}

/*void LocalAlignment::read_config(const string &fileName)
{
	ifstream cFile(fileName);
	if (cFile.is_open())
	{
		string line;
		while (getline(cFile, line)) {
			line.erase(remove_if(line.begin(), line.end(), isspace),
				line.end());
			if (line[0] == '#' || line.empty())
				continue;
			unsigned int delimiterPos = line.find("=");
			string name = line.substr(0, delimiterPos);
			int value = stoi(line.substr(delimiterPos + 1));
			if (name == "MATCH") this->match = value;
			if (name == "MISMATCH") this->mismatch = value;
			if (name == "GAP") this->gap = value;
			if (name == "THRESHOLD") this->threshold = value;
		}

	}
	else {
		cerr << "Couldn't open config file for reading.\n";
	}
}*/
