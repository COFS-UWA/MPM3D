#ifndef __Stack_List_h__
#define __Stack_List_h__

// Stack like list
/*
Explanation of offset:
For any class Item, add member Item *next (can be other name)
	class Item
	{
	....
	Item *next;
	}

offset is offsetof(Item, next)
*/
template <typename Item, size_t offset>
class StackList
{
protected:
	Item *top;
	size_t num;

	inline static Item *get_item_from_pointer(Item **pt)
	{
		return (Item *)((char *)pt - offset);
	}
	inline static Item **get_pointer_from_item(Item *it)
	{
		return (Item **)((char *)it + offset);
	}

public:
	StackList() : top(nullptr), num(0) {}
	~StackList() { reset(); }

	inline void push(Item *it)
	{
		*get_pointer_from_item(it) = top;
		top = it;
		++num;
	}

	inline void push(Item &it) { push(&it); }

	inline Item *pop()
	{
		Item *last_item = top;
		Item *sec_last_item = *get_pointer_from_item(top);
		top = sec_last_item;
		--num;
		return last_item;
	}
	
	inline void reset()
	{
		top = nullptr;
		num = 0;
	}

	inline size_t get_num() { return num; }
	inline bool is_empty() { return num == 0; }

	inline Item *first() { return top; }
	inline static Item *next(Item *it) { return *get_pointer_from_item(it); }
	inline bool static is_end(Item *it) { return it == nullptr; }
	inline bool static is_not_end(Item *it) { return it != nullptr; }
};

#endif