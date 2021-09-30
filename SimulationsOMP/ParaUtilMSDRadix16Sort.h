#ifndef __Para_Util_MSD_Radix_16_Sort_h__
#define __Para_Util_MSD_Radix_16_Sort_h__

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#include <tbb/tbb.h>

namespace ParaUtil
{
namespace Internal
{
    constexpr size_t msd_radix_16_bit_num = 4;
    constexpr size_t min_data_num_for_msd_radix_sort_16 = 1024;

    // initial key_bit_num shouldn't be 0
    class MSDRadixSort16Task : public tbb::task
    {
    public:
        MSDRadixSort16Task(size_t _keys[], size_t _values[],
            size_t _data_num, size_t _key_bit_num,
            size_t _tmp_keys[], size_t _tmp_values[]) :
            keys(_keys), values(_values), data_num(_data_num),
            key_bit_num(_key_bit_num),
            tmp_keys(_tmp_keys), tmp_values(_tmp_values),
            is_continuation(false) { }
        tbb::task* execute();

    protected:
        size_t *keys, *values;
        size_t data_num, key_bit_num;
        size_t *tmp_keys, *tmp_values;
        bool is_continuation;
    };

    void msd_radix_sort_16(size_t keys[], size_t values[],
        size_t data_num, size_t max_key,
        size_t tmp_keys[], size_t tmp_values[],
        int thread_num);
}
}

#endif