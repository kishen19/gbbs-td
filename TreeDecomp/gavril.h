#include "gbbs/gbbs.h"

namespace gbbs{

template <class Graph>
size_t GavrilTreeDecomp(Graph& GA, sequence<uintE>& pi, sequence<sequence<uintE>>& B, sequence<uintE>& T, sequence<size_t>& bag_size){
    using W = typename Graph::weight_type;
	size_t n = GA.n;

	auto pi_inv = sequence<uintE>::uninitialized(n);
    parallel_for(0, n, [&](size_t i) {pi_inv[pi[i]] = i;});

    T[n] = 0;
    parallel_for(0, n, [&](size_t i) {
        uintE pos_u = pi_inv[i];
        auto count_f = [&](const uintE& u, const uintE& v, const W& wgh) {
            uintE pos_v = pi_inv[v];
            return pos_u < pos_v;
        };
        bag_size[pos_u] = GA.get_vertex(i).out_neighbors().count(count_f);
        auto filter_f = [&](const uintE& u, const uintE& v) {
            uintE pos_v = pi_inv[v];
            return pos_u < pos_v;
        };
        B[pos_u] = GA.get_vertex(i).out_neighbors().filter2(filter_f);
        auto sort_f = [&](uintE v){return pi_inv[v];};
        parlay::integer_sort_inplace(B[pos_u],sort_f);
        T[pos_u] = B[pos_u][0];
    });
    return parlay::reduce_max(bag_size);
}

}