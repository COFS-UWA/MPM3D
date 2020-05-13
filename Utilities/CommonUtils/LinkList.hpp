#ifndef __Link_List_h__
#define __Link_List_h__

// Double link list
// LinkListPointer should be member of class Item
// offset: offsetof(Item, pointer_member)
template <typename Item>
struct LinkListPointer
{
	Item *prev;
	Item *next;
};

template <typename Item, size_t offset>
class LinkList
{
protected:
	typedef LinkListPointer<Item> Pointer;
	
	size_t num;
	const Item *head_item;
	Pointer head_pointer;

	inline static Item *get_item_from_pointer(Pointer *pointer)
	{
		return (Item *)((char *)pointer - offset);
	}
	inline static Pointer *get_pointer_from_item(Item *item)
	{
		return (Pointer *)((char *)item + offset);
	}

public:
	LinkList() : num(0),
		head_item(get_item_from_pointer(&head_pointer))
	{
		head_pointer.prev = const_cast<Item *>(head_item);
		head_pointer.next = const_cast<Item *>(head_item);
	}
	~LinkList() { reset(); }

	inline void append(Item *item)
	{
		Item *last_item = head_pointer.prev;
		Pointer *pointer = get_pointer_from_item(item);
		pointer->prev = last_item;
		pointer->next = const_cast<Item *>(head_item);
		get_pointer_from_item(last_item)->next = item;
		head_pointer.prev = item;
		++num;
	}

	inline void prepend(Item *item)
	{
		Item *first_item = head_pointer.next;
		Pointer *pointer = get_pointer_from_item(item);
		pointer->prev = const_cast<Item *>(head_item);
		pointer->next = first_item;
		get_pointer_from_item(first_item)->prev = item;
		head_pointer.next = item;
		++num;
	}

	inline void del(Item *item)
	{
		Pointer *pointer = get_pointer_from_item(item);
		get_pointer_from_item(pointer->next)->prev = pointer->prev;
		get_pointer_from_item(pointer->prev)->next = pointer->next;
		--num;
	}

	inline void append(Item &item) { append(&item); }
	inline void prepend(Item &item) { prepend(&item); }
	inline void del(Item &item) { del(&item); }

	// push the item into the first place
	inline void push(Item *item) { prepend(item); }
	inline void push(Item &item) { push(&item); }
	
	// pop the first item
	inline Item *pop()
	{
		Item *res = first();
		del(res);
		return res;
	}

	// insert new_item before in_list_item
	inline void insert_before(Item *in_list_item, Item *new_item)
	{
		Pointer *in_list_pt = get_pointer_from_item(in_list_item);
		Item *prev_item = in_list_pt->prev;
		Pointer *new_pointer = get_pointer_from_item(new_item);
		new_pointer->prev = prev_item;
		new_pointer->next = in_list_item;
		get_pointer_from_item(prev_item)->next = new_item;
		in_list_pt-> prev = new_item;
		++num;
	}

	// insert new_item after in_list_item
	inline void insert_after(Item *in_list_item, Item *new_item)
	{
		Pointer *in_list_pt = get_pointer_from_item(in_list_item);
		Item *next_item = in_list_pt->next;
		Pointer *new_pointer = get_pointer_from_item(new_item);
		new_pointer->prev = in_list_item;
		new_pointer->next = next_item;
		in_list_pt->next = new_item;
		get_pointer_from_item(next_item)->prev = new_item;
		++num;
	}

	inline void insert_before(Item &in_list_item, Item &new_item)
	{
		insert_before(&in_list_item, &new_item);
	}
	inline void insert_after(Item &in_list_item, Item &new_item)
	{
		insert_after(&in_list_item, &new_item);
	}

	inline void reset()
	{
		head_pointer.prev = const_cast<Item *>(head_item);
		head_pointer.next = const_cast<Item *>(head_item);
		num = 0;
	}

	// transfer ownership from another link to this list
	inline void transfer(LinkList &another)
	{
		if (another.is_empty())
		{
			reset();
			return;
		}
		Pointer &an_hpt = another.head_pointer;
		head_pointer.prev = an_hpt.prev;
		head_pointer.next = an_hpt.next;
		get_pointer_from_item(an_hpt.prev)->next = const_cast<Item *>(head_item);
		get_pointer_from_item(an_hpt.next)->prev = const_cast<Item *>(head_item);
		num = another.num;
		another.reset();
	}

	inline size_t get_num() { return num; }
	inline bool is_empty() { return num == 0; }

	inline Item *first() { return head_pointer.next; }
	inline Item *last() { return head_pointer.prev; }
	inline static Item *next(Item *item) { return get_pointer_from_item(item)->next; }
	inline static Item *prev(Item *item) { return get_pointer_from_item(item)->prev; }
	inline bool is_end(Item *item) { return item == const_cast<Item *>(head_item); }
	inline bool is_not_end(Item *item) { return item != const_cast<Item *>(head_item); }
};

#endif