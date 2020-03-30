#pragma once

struct Coordinate
{
	unsigned int row;
	unsigned int col;
	inline Coordinate()
	{

	}
	inline Coordinate(const unsigned int &row, const unsigned int &col)
	{
		this->row = row;
		this->col = col;
	}
	inline Coordinate &operator= (const Coordinate &c)
	{
		this->row = c.row;
		this->col = c.col;
		return *this;
	}
	inline bool operator== (const Coordinate &c) const
	{
		return this->row == c.row && this->col == c.col;
	}
	
	inline bool operator!= (const Coordinate &c) const
	{
		return !(this->row == c.row && this->col == c.col);
	}
	inline bool operator< (const Coordinate &c) const
	{
		if (this->col == c.col) return this->row < c.row;
		return this->col < c.col;
	}
	inline bool isNull() const
	{
		return this->row == 0 || this->col == 0;
	}
};

namespace std
{
	template <> struct std::hash<Coordinate> {
		inline size_t operator()(const Coordinate &v) const {
			std::hash<int> int_hasher;
			return int_hasher(v.row) ^ int_hasher(v.col);
		}
	};
};

struct Coordinate_Score
{
	Coordinate c;
	unsigned int score;
	inline Coordinate_Score()
	{

	}
	inline Coordinate_Score(const Coordinate &c, const unsigned int &score)
	{
		this->c = c;
		this->score = score;
	}
	inline bool operator==(const Coordinate_Score &cs) const
	{
		return this->c == cs.c && this->score == cs.score;
	}
	inline bool operator<(const Coordinate_Score &cs) const
	{
		return this->score < cs.score;
	}
	inline bool operator>(const Coordinate_Score &cs) const
	{
		return this->score > cs.score;
	}
};

struct Interval_Coordinate
{
	Coordinate start;
	Coordinate end;
	unsigned int score;
	inline Interval_Coordinate()
	{

	}
	inline Interval_Coordinate(const Coordinate &start, const Coordinate &end, const unsigned int &score)
	{
		this->start = start;
		this->end = end;
		this->score = score;
	};
	inline bool operator==(const Interval_Coordinate& ic) const
	{
		return this->start == ic.start && this->end == ic.end && this->score == ic.score;
	}
};



