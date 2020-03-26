#ifndef __Memory_Utils_String_Buffer_h__
#define __Memory_Utils_String_Buffer_h__

#include <string>

namespace MemoryUtils
{
	/*=======================================================
	 * Class StringBuffer
	 *-------------------------------------------------------
	 *   A buffer for C styel string. 
	 * Note: 
	 *   1. The member variable "length" doesn't include '\0'
	 *      at the end, but member variable "capacity" take
	 *      it into consideration.
	 *======================================================*/
	class StringBuffer
	{
	protected:
		char *str;
		size_t length;
		size_t capacity;
	public:
		StringBuffer(size_t init_string_length = 0) :
			str(nullptr), length(0), capacity(0)
		{
			if (init_string_length && (++init_string_length) > capacity)
				reserve_raw(init_string_length);
		}
		~StringBuffer() { if (str) delete[] str; }

		inline const char *get_string(void) const { return str; }
		inline const char *c_str(void) const { return str; }
		inline size_t get_length(void) const { return length; }
		inline size_t get_capacity(void) const { return capacity; }

		inline void reserve(size_t reserve_string_length)
		{
			if (reserve_string_length && (++reserve_string_length) > capacity)
				reserve_raw(reserve_string_length);
		}
		inline void expand(size_t expand_string_length)
		{
			if (expand_string_length)
				reserve_raw(capacity + expand_string_length);
		}
		inline void reset(void)
		{
			length = 0;
			if (str) str[0] = '\0';
		}
		void clear(void)
		{
			if (str)
			{
				delete[] str;
				str = nullptr;
			}
			length = 0;
			capacity = 0;
		}

		void add_string(const char *another_str)
		{
			size_t new_capacity = length + strlen(another_str) + 1;
			if (new_capacity > capacity)
				reserve_raw(new_capacity);
			strcpy(str + length, another_str);
			length = new_capacity - 1;
		}
		void add_char(char another_char)
		{
			size_t new_capacity = length + 2;
			if (new_capacity > capacity)
				reserve_raw(new_capacity);
			str[length] = another_char;
			++length;
			str[length] = '\0';
		}

	protected:
		void reserve_raw(size_t reserve_size)
		{
			char *tmp = new char[reserve_size];
			tmp[0] = '\0';
			if (str)
			{
				if (length)
					strcpy(tmp, str);
				delete[] str;
			}
			str = tmp;
			capacity = reserve_size;
		}

	public:
		inline StringBuffer &operator<< (StringBuffer &another_buf)
		{
			add_string(another_buf.get_string());
			return *this;
		}
		inline StringBuffer &operator<< (const char *another_str)
		{
			add_string(another_str);
			return *this;
		}
		inline StringBuffer &operator<< (char another_char)
		{
			add_char(another_char);
			return *this;
		}
		inline StringBuffer &operator<< (std::string &another_str)
		{
			add_string(another_str.c_str());
			return *this;
		}
		inline StringBuffer &operator<< (size_t num)
		{
			char str_tmp[25];
			snprintf(str_tmp, 25, "%llu", num);
			add_string(str_tmp);
			return *this;
		}
		inline bool operator== (StringBuffer &another_buf)
		{
			return strcmp(str, another_buf.str) == 0;
		}
		inline bool operator!= (StringBuffer &another_buf)
		{
			return strcmp(str, another_buf.str) != 0;
		}
		inline bool operator== (const char *another_str)
		{
			return strcmp(str, another_str) == 0;
		}
		inline bool operator!= (const char *another_str)
		{
			return strcmp(str, another_str) != 0;
		}
		inline bool operator== (std::string &another_str)
		{
			return strcmp(str, another_str.c_str()) == 0;
		}
		inline bool operator!= (std::string &another_str)
		{
			return strcmp(str, another_str.c_str()) != 0;
		}
	};

	std::ostream &operator<< (std::ostream& out, StringBuffer &str_buf)
	{
		out << str_buf.get_string();
		return out;
	}
};

#endif