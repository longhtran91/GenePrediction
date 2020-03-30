#include "ExonChaining.h"

using namespace std;

ExonChaining::ExonChaining()
{

}

ExonChaining::ExonChaining(const list<Interval_Coordinate> &exons)
{
	//create ordered vertices vector & create right2leftweight multimap
	multimap<Coordinate, Coordinate_Score> rightC2leftCweight;
	multimap<unsigned int, Coordinate_Score> right2leftweight;	
	vector<Coordinate> ordered_vertices;
	for (const auto &elem : exons)
	{
		rightC2leftCweight.insert({ elem.end, Coordinate_Score(elem.start, elem.score) });
		right2leftweight.insert({ elem.end.col, Coordinate_Score(elem.start, elem.score) });
		ordered_vertices.push_back(elem.start);
		ordered_vertices.push_back(elem.end);
	}
	sort(ordered_vertices.begin(), ordered_vertices.end());
	ordered_vertices.erase(unique(ordered_vertices.begin(), ordered_vertices.end()), ordered_vertices.end());
	unsigned int size = ordered_vertices.size();
	ordered_vertices.reserve(size);
	vector<unsigned int> scores(size);
	
	unsigned int pre_j = 0;
	for (unsigned int i = 1; i < size; ++i)
	{
		unsigned int right_j = ordered_vertices[i].col;
		if (right_j == pre_j)
		{ 
			scores[i] = scores[i - 1];
			continue;
		}

		const auto &itr_right = right2leftweight.equal_range(right_j);
		if (itr_right.first == itr_right.second) //current index is left end
		{
			scores[i] = scores[i - 1];
		}
		else //current index is right
		{
			//find all possible lefts
			auto itr_left = itr_right.first;
			unsigned int max_weight = 0;
			while (itr_left != itr_right.second)
			{
				unsigned int left = itr_left->second.c.col;
				unsigned int itv_weight = itr_left->second.score;
				const auto &left_itr_index = lower_bound(ordered_vertices.cbegin(), ordered_vertices.cbegin() + i, Coordinate(0, left),
																											[](const Coordinate &lhs, Coordinate const & rhs)
																											{return lhs.col < rhs.col; });
				unsigned int left_index = distance(ordered_vertices.cbegin(), left_itr_index);
				if (scores[left_index] + itv_weight > max_weight)
				{
					max_weight = scores[left_index] + itv_weight;
				}
				++itr_left;
			}
			scores[i] = max({ max_weight, scores[i - 1] });
		}
		pre_j = right_j;
	}
	
	//backtrack scores to get paths
	int index = size - 1;
	Coordinate right, left;
	unsigned int weight;
	while (index >= 0)
	{
		right = ordered_vertices[index];
		const auto &itr_right = rightC2leftCweight.find(right); //find all the lefts-weights of current right
		if (itr_right == rightC2leftCweight.cend()) //if the current index is left
		{
			//look for any right end at the same j
			const auto &itr_right = right2leftweight.find(right.row);
			if (itr_right == right2leftweight.cend()) //no right end found
			{
				--index;
			}
			else //right end at the same j found
			{
				unsigned int tempw = scores[index];			
				
				//find all coordinates including left & right at this j
				auto itr_cor = equal_range(ordered_vertices.cbegin(), ordered_vertices.cend(), Coordinate(0, right.row), [](const Coordinate &lhs, Coordinate const & rhs)
																													{return lhs.col < rhs.col; });
				Coordinate new_right;
				Coordinate new_left;
				unsigned int max_itv_weight = 0;
				auto cor = itr_cor.first;
				while (itr_cor.first != itr_cor.second)
				{
					auto itr_rights = rightC2leftCweight.equal_range(*itr_cor.first); //find all the rights
					while (itr_rights.first != itr_rights.second) //find rights with highest interval score
					{
						if (itr_rights.first->second.score > max_itv_weight)
						{
							cor = itr_cor.first;
							new_right = itr_rights.first->first;
							new_left = itr_rights.first->second.c;
							max_itv_weight = itr_rights.first->second.score;
						}
						++itr_rights.first;
					}
					++itr_cor.first;

				}

				//const auto &itr_new_right = rightC2leftCweight.find(new_right);
				unsigned int new_weight = max_itv_weight;
                //const auto &new_left_itr = lower_bound(ordered_vertices.cbegin(), cor, new_left);
				unsigned int new_left_index = distance(ordered_vertices.cbegin(), cor);
				if (scores[new_left_index] + new_weight == tempw)
				{
					//add the interval
					this->intervals.push_front(Interval_Coordinate(new_left, new_right, new_weight));
					index = new_left_index; //go to left as next
				}
				else
				{
					--index;;
				}
			}
				
		}
		else //if the index is right
		{
			//find max score of this right
			left = itr_right->second.c;
			weight = itr_right->second.score;
			const auto &left_itr = lower_bound(ordered_vertices.cbegin(), ordered_vertices.cbegin() + index, left);

			unsigned int left_index = distance(ordered_vertices.cbegin(), left_itr);
			unsigned int test_score = scores[left_index] + weight; //test_score to test it against the score in the scores vector
			unsigned int score = scores[index]; //score in the scores vector
			if (test_score == score) // if test_score = score at this position, use this interval
			{
				//add the interval
				this->intervals.push_front(Interval_Coordinate(left, right, weight));
				index = left_index; //go to left as next
			}
			else
			{
				--index;
			}
		}
	}
}

ExonChaining::~ExonChaining()
{

}

list<Interval_Coordinate> ExonChaining::get_intervals()
{
	return this->intervals;
}
