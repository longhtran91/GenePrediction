#pragma once
#include <array>
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include <sstream>
#include <iostream>

#define DIAGONAL_PATH 'd'
#define VERTICAL_PATH 'u'
#define HORIZONTAL_PATH 'l'

using namespace std;
namespace std
{
	template<typename T, size_t N>
	struct hash<array<T, N> >
	{
		typedef array<T, N> argument_type;
		typedef size_t result_type;

		result_type operator()(const argument_type& a) const
		{
			hash<T> hasher;
			result_type h = 0;
			for (result_type i = 0; i < N; ++i)
			{
				h = h * 31 + hasher(a[i]);
			}
			return h;
		}
	};
}

class LocalAlignment
{

private:

	unordered_set<array<unsigned, 3>> exons;
	string dna_sequence;
	string dna_template;
	unsigned int i_size, j_size, threshold;
	int match, mismatch, gap;

	string read_file(const string &fileName);
	void read_config(const string &fileName);


public:
	LocalAlignment();
	LocalAlignment(const string &dna_sequence_filename, const string &dna_template_filename);
	~LocalAlignment();
	unordered_set<array<unsigned, 3>> getExons() const;
};
