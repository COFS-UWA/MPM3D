#ifndef __Num_Pair_Hash_Table_hpp__
#define __Num_Pair_Hash_Table_hpp__

#include "ItemBuffer.hpp"

// assistent data structs for tetrahedron mesh
template <typename Num>
class NumPairHashTable
{
public:
	struct NumPair
	{
		Num key;
		Num val;
		NumPair *next;
	};

protected:
	size_t entry_num;
	NumPair **num_pair_list;
	MemoryUtils::ItemBuffer<NumPair> num_pair_buf;
	size_t pair_num;

public:
	NumPairHashTable(size_t _entry_num) : 
		entry_num(_entry_num), pair_num(0)
	{
		num_pair_list = new NumPair*[entry_num];
		for (size_t en_id = 0; en_id < entry_num; ++en_id)
			num_pair_list[en_id] = nullptr;
	}
	~NumPairHashTable() {}

	size_t get_pair_num() { return pair_num; }

	template <typename ResNum = Num>
	void output_pairs(ResNum *mem)
	{
		for (size_t e_id = 0; e_id < entry_num; ++e_id)
			for (NumPair *pair = num_pair_list[e_id]; pair; pair = pair->next)
			{
				mem[0] = ResNum(pair->key);
				mem[1] = ResNum(pair->val);
				mem += 2;
			}
	}

	// Assume key has single value
	bool find_key(Num _key, Num &_val)
	{
		for (NumPair *pair = num_pair_list[_key % entry_num];
			 pair; pair = pair->next)
			if (pair->key == _key)
			{
				_val = pair->val;
				return true;
			}
		return false;
	}
	
	bool add_key(Num _key, Num &_val)
	{
		NumPair *&top = num_pair_list[_key % entry_num];
		for (NumPair *pair = top; pair; pair = pair->next)
		{
			if (pair->key == _key)
				return false;
		}
		NumPair *tmp = num_pair_buf.alloc();
		tmp->key = _key;
		tmp->val = _val;
		tmp->next = top;
		top = tmp;
		++pair_num;
		return true;
	}

	// Assume key has multiple value
	bool find_pair(Num _key, Num _val)
	{
		for (NumPair *pair = num_pair_list[_key % entry_num];
			 pair; pair = pair->next)
			if (pair->key == _key && pair->val == _val)
				return true;
		return false;
	}

	// return false if the pair already exists
	bool add_pair(Num _key, Num _val)
	{
		NumPair *&top = num_pair_list[_key % entry_num];
		for (NumPair *pair = top; pair; pair = pair->next)
		{
			if (pair->key == _key && pair->val == _val)
				return false;
		}
		NumPair *tmp = num_pair_buf.alloc();
		tmp->key = _key;
		tmp->val = _val;
		tmp->next = top;
		top = tmp;
		++pair_num;
		return true;
	}

	void reset()
	{
		for (size_t en_id = 0; en_id < entry_num; ++en_id)
			num_pair_list[en_id] = nullptr;
		num_pair_buf.reset();
		pair_num = 0;
	}
};

#endif