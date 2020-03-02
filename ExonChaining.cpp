#include "ExonChaining.h"

using namespace std;

ExonChaining::ExonChaining()
{

}

ExonChaining::ExonChaining(const unordered_set<array<unsigned int, 3>> &exons)
{

	multimap<unsigned int, array<unsigned int, 2>> right2leftweight;

	//create ordered vertices vector & create right2leftweight multimap
	set<unsigned int> temp_set; // to temporary store all the vertices to a set
	for (const auto &elem : exons) {
		right2leftweight.insert(pair<unsigned int, array<unsigned int, 2>>(elem[1], { elem[0],  elem[2] }));
		temp_set.insert(elem[0]);
		temp_set.insert(elem[1]);
	}
	unsigned int size = temp_set.size();

	vector<unsigned int> ordered_vertices;
	ordered_vertices.reserve(size);
	ordered_vertices.assign(temp_set.begin(), temp_set.end());
	unsigned int *scores = new unsigned int[size];  //create ordered verticies vector from set to remove duplicate and sort by default


	//fill scores vector
	for (unsigned int i = 1; i < size; i++)
	{
		unsigned int right = (ordered_vertices)[i];

		//find highest score of multiples intervals with the same rights
		pair <multimap<unsigned int, array<unsigned int, 2>>::const_iterator, multimap<unsigned int, array<unsigned int, 2>>::const_iterator> itr_right2lw; //pair of multimap iterators
		itr_right2lw = right2leftweight.equal_range(right); //find all the lefts-weights of current right

		if (itr_right2lw.first == itr_right2lw.second) //if the index is left
		{
			(scores)[i] = (scores)[i - 1];
		}
		else //if the index is right
		{
			multimap<unsigned int, array<unsigned int, 2>>::const_iterator max_itr = find_max_itr(itr_right2lw.first, itr_right2lw.second);
			unsigned int left = max_itr->second[0];
			unsigned int weight = max_itr->second[1];

			//find the left index
			unsigned int left_index = lower_bound((ordered_vertices).begin(), (ordered_vertices).end(), left) - (ordered_vertices).begin();

			//find max score of all possible scores
			unsigned int possible_scores[2] = { (scores)[left_index] + weight , (scores)[i - 1] };
			(scores)[i] = *max_element(possible_scores, possible_scores + 2);
		}
	}

	//backtrack scores to get paths
	int index = size - 1;
	unsigned int right, left, weight;

	while (index >= 0)
	{
		right = (ordered_vertices)[index];
		pair <multimap<unsigned int, array<unsigned int, 2>>::const_iterator, multimap<unsigned int, array<unsigned int, 2>>::const_iterator> itr_right2lw; //pair of multimap iterators
		itr_right2lw = right2leftweight.equal_range(right); //find all the lefts-weights of current right
		if (itr_right2lw.first == itr_right2lw.second) //if the index is left
		{
			index--;
		}
		else //if the index is right
		{
			multimap<unsigned int, array<unsigned int, 2>>::const_iterator max_itr = find_max_itr(itr_right2lw.first, itr_right2lw.second);
			left = max_itr->second[0];
			weight = max_itr->second[1];
			unsigned left_index = lower_bound((ordered_vertices).begin(), (ordered_vertices).end(), left) - (ordered_vertices).begin();
			unsigned int test_score = (scores)[left_index] + weight; //test_score to test it against the score in the scores vector
			unsigned int score = (scores)[index]; //score in the scores vector
			if (test_score == score) // if test_score = score at this position, use this interval
			{
				//add the interval
				this->intervals.push({ left, right, weight });
				index = left_index; //go to left as next
			}
			else
			{
				index--;
			}
		}
	}
    delete[] scores;
}

ExonChaining::~ExonChaining()
{

}

template <typename T>
T ExonChaining::find_max_itr(T start, T end)
{
	T max_itr = start;
	unsigned int max_score = max_itr->second[1];
	for (T itr = start; itr != end; itr++)
	{
		if (itr->second[1] > max_score)
		{
			max_score = itr->second[1];
			max_itr = itr;
		}
	}
	return max_itr;
}

stack<array<unsigned int, 3>> ExonChaining::get_intervals()
{
	return this->intervals;
}
