#include "Tests_pcp.h"

#include "LinkList.hpp"
#include "StackList.hpp"

namespace
{
	struct IntItem
	{
		int num;
		LinkListPointer<IntItem> pointer;
		char padding1;
		IntItem *next_by_stack;
	};

	std::ostream &operator<<(std::ostream &out, IntItem &it)
	{
		std::cout << it.num;
		return out;
	}

	template <typename Item, size_t offset>
	void print_list(LinkList<Item, offset> &llist)
	{
		for (Item *it_iter = llist.first();
			llist.is_not_end(it_iter);
			it_iter = llist.next(it_iter))
			std::cout << size_t(it_iter)
					  << " " << size_t(it_iter->pointer.prev)
					  << " " << size_t(it_iter->pointer.next)
					  << " " << *it_iter
					  << " num: " << llist.get_num() << "\n";
	}

	template <typename Item, size_t offset>
	void print_list(StackList<Item, offset> &slist)
	{
		for (Item *it_iter = slist.first(); 
			slist.is_not_end(it_iter);
			it_iter = slist.next(it_iter))
			std::cout << size_t(it_iter) << " "
					  << size_t(it_iter->next_by_stack) << " "
					  << *it_iter
					  << " num: " << slist.get_num() << "\n";
	}
}

void test_stack_and_link_list()
{
	std::cout << "test link list:\n";

	LinkList<IntItem, offsetof(IntItem, pointer)> llist;

	std::cout << "test append\n";
	IntItem items[10];
	for (size_t i = 0; i < 10; i++)
	{
		items[i].num = i;
		llist.append(items[i]);
	}

	//std::cout << offsetof(IntItem, pointer) << " ";
	//LinkListPointer<IntItem> *pt = llist.get_pointer_from_item(items);
	//std::cout << size_t(items) << " " 
	//		  << size_t(pt) << " " 
	//		  << size_t(&(items->pointer)) << " "
	//		  << size_t(llist.get_item_from_pointer(pt)) << "\n";

	print_list(llist);

	std::cout << "test reset\n";
	llist.reset();
	if (llist.is_end(llist.first()) && llist.is_empty())
		std::cout << "list is empty after being reset.\n";
	
	std::cout << "test prepend\n";
	for (size_t i = 0; i < 10; i++)
	{
		items[i].num = i + 1;
		llist.prepend(items[i]);
	}

	print_list(llist);
	
	std::cout << "test del\n";
	llist.del(items[0]);
	llist.del(items[2]);
	llist.del(items[9]);
	llist.del(items[8]);

	print_list(llist);

	std::cout << "test stack list\n";
	StackList<IntItem, offsetof(IntItem, next_by_stack)> slist;

	std::cout << "test push\n";
	for (size_t i = 0; i < 10; i++)
		slist.push(items[i]);

	print_list(slist);

	std::cout << "test pop\n";
	for (size_t i = 0; i < 10; i++)
	{
		IntItem *it = slist.pop();
		std::cout << *it << " num: " << slist.get_num() << "\n";
	}

	std::cout << "test reset\n";
	slist.push(items[2]);
	slist.push(items[1]);
	slist.reset();
	if (slist.is_end(slist.first()) && slist.is_empty())
		std::cout << "list is empty after being reset.\n";

}
