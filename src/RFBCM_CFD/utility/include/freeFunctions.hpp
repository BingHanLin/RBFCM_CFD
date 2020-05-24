#include <algorithm>  // std::sort, std::stable_sort
#include <numeric>    // std::iota
#include <vector>

using namespace std;

template <typename T>
std::vector<size_t> sortIndexes(const std::vector<T> &vec)
{
    // initialize original index locations
    std::vector<size_t> idx(vec.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    std::stable_sort(idx.begin(), idx.end(), [&vec](size_t i1, size_t i2) {
        return vec[i1] < vec[i2];
    });

    return idx;
}