#include "LocalAlignment.h"

using namespace std;

LocalAlignment::LocalAlignment()
{

}

LocalAlignment::LocalAlignment(const string &dna_sequence_filename, const string &dna_template_filename)
{
	this->read_config("var.config");
	this->dna_template = read_file(dna_template_filename);
	this->dna_sequence = read_file(dna_sequence_filename);

	this->i_size = dna_template.length() + 1;
	this->j_size = dna_sequence.length() + 1;
	int **scores_matrix = new int *[this->i_size];
	unordered_set<array<unsigned int, 2>> **backtracking_matrix = new unordered_set<array<unsigned int, 2>> *[this->i_size];
	
	//fill the rest of matrix by dynamic programming local alignment
	for (unsigned int i = 0; i < this->i_size; i++)
	{
		scores_matrix[i] = new int[this->j_size];
		backtracking_matrix[i] = new unordered_set<array<unsigned int, 2>>[this->j_size];

		//rolling delete used/obsolete rows
		if (i >= 2)
		{
				delete[] scores_matrix[i-2];
				delete[] backtracking_matrix[i - 2];
		}

		for (unsigned int j = 0; j < this->j_size; j++)
		{
			//initialize first row and column to be 0
			if (i == 0 || j == 0)
			{
				scores_matrix[i][j] = 0;
				backtracking_matrix[i][j].insert({ i, j });
				continue; 
			}

			bool free_ride = true;
			int vertical_score, horizontal_score, diagonal_score;

			horizontal_score = scores_matrix[i][j-1] + this->gap;
			vertical_score = scores_matrix[i-1][j] + this->gap;

			if (dna_sequence[j - 1] == dna_template[i - 1])	diagonal_score = scores_matrix[i - 1][j - 1] + this->match;
			else diagonal_score = scores_matrix[i - 1][j - 1] + this->mismatch;

            int scores[4] = {0, vertical_score, horizontal_score, diagonal_score };
			unsigned int score = *max_element(scores, scores + 4);

			if (score == vertical_score)
			{
				free_ride = false;
				backtracking_matrix[i][j].insert(backtracking_matrix[i-1][j].cbegin(), backtracking_matrix[i - 1][j].cend());
			}

			if (score == horizontal_score)
			{
				free_ride = false;
				backtracking_matrix[i][j].insert(backtracking_matrix[i][j-1].cbegin(), backtracking_matrix[i][j-1].cend());
			}

			if (score == diagonal_score)
			{
				free_ride = false;
				backtracking_matrix[i][j].insert(backtracking_matrix[i-1][j - 1].cbegin(), backtracking_matrix[i-1][j - 1].cend());
			}

			if (free_ride)  backtracking_matrix[i][j].insert({ i, j });

			scores_matrix[i][j] = score;
			
			if (score >= this->threshold)
			{
				for (const auto& elem : backtracking_matrix[i][j]) {
					this->exons.insert({ (elem[1] == 0) ? 1 : elem[1],j, score});
				}
			}
		}
	}

	for (int k = 1; k < 3; k++)
	{
		delete[] scores_matrix[this->i_size - k];
		delete[] backtracking_matrix[this->i_size - k];
	}
	delete[] scores_matrix;
	delete[] backtracking_matrix;
}

LocalAlignment::~LocalAlignment()
{

}

unordered_set<array<unsigned, 3>> LocalAlignment::getExons() const
{
	return this->exons;
}

void LocalAlignment::read_config(const string &fileName)
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
}

string LocalAlignment::read_file(const string &fileName)
{
	string data;
	string line;
	ifstream ifs(fileName.c_str());

	while (getline(ifs, line))
	{
		if (line[0] == '>') continue;
		data += line;
	}
	return data;
}
