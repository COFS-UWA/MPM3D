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

	inline static Item *get_item_from_pointer(Pointer *pt)
	{
		return (Item *)((char *)pt - offset);
	}
	inline static Pointer *get_pointer_from_item(Item *it)
	{
		return (Pointer *)((char *)it + offset);
	}

public:
	LinkList() : num(0),
		head_item(get_item_from_pointer(&head_pointer))
	{
		head_pointer.prev = const_cast<Item *>(head_item);
		head_pointer.next = const_cast<Item *>(head_item);
	}
	~LinkList() { reset(); }

	inline void append(Item *it)
	{
		Item *last_item = head_pointer.prev;
		Pointer *pt = get_pointer_from_item(it);
		pt->prev = last_item;
		pt->next = const_cast<Item *>(head_item);
		get_pointer_from_item(last_item)->next = it;
		head_pointer.prev = it;
		++num;
	}

	inline void prepend(Item *it)
	{
		Item *first_item = head_pointer.next;
		Pointer *pt = get_pointer_from_item(it);
		pt->prev = const_cast<Item *>(head_item);
		pt->next = first_item;
		get_pointer_from_item(first_item)->prev = it;
		head_pointer.next = it;
		++num;
	}

	inline void del(Item *it)
	{
		Pointer *pt = get_pointer_from_item(it);
		get_pointer_from_item(pt->next)->prev = pt->prev;
		get_pointer_from_item(pt->prev)->next = pt->next;
		--num;
	}

	inline void append(Item &it) { append(&it); }
	inline void prepend(Item &it) { prepend(&it); }
	inline void del(Item &it) { del(&it); }

	inline void reset()
	{
		head_pointer.prev = const_cast<Item *>(head_item);
		head_pointer.next = const_cast<Item *>(head_item);
		num = 0;
	}

	inline size_t get_num() { return num; }
	inline bool is_empty() { return num == 0; }

	inline Item *first() { return head_pointer.next; }
	inline Item *last() { return head_pointer.prev; }
	inline static Item *next(Item *it) { return get_pointer_from_item(it)->next; }
	inline static Item *prev(Item *it) { return get_pointer_from_item(it)->prev; }
	inline bool is_end(Item *it) { return it == const_cast<Item *>(head_item); }
	inline bool is_not_end(Item *it) { return it != const_cast<Item *>(head_item); }
};

#endif