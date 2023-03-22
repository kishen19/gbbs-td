#pragma once

#include "gbbs/gbbs.h"

namespace gbbs{

template <class Graph>
size_t GavrilTreeDecomp(Graph& GA, sequence<uintE>& pi, sequence<sequence<uintE>>& B, sequence<uintE>& T, sequence<size_t>& bag_size){
    using W = typename Graph::weight_type;
	size_t n = GA.n;

	auto pi_inv = sequence<uintE>::uninitialized(n);
    parallel_for(0, n, [&](size_t i) {pi_inv[pi[i]] = i;});

    T[n] = n;
    parallel_for(0, n, [&](uintE i) { // vertex i
        uintE pos_u = pi_inv[i];
        auto count_f = [&](const uintE& u, const uintE& v, const W& wgh) {
            uintE pos_v = pi_inv[v];
            return pos_u < pos_v;
        };
        bag_size[i] = GA.get_vertex(i).out_neighbors().count(count_f);
        if (bag_size[i] == 0){
            T[i] = n;
        } else{
            auto filter_f = [&](const uintE& u, const uintE& v) {
                uintE pos_v = pi_inv[v];
                return pos_u < pos_v;
            };
            B[i] = GA.get_vertex(i).out_neighbors().filter2(filter_f, bag_size[i]);
            auto sort_f = [&](uintE v){return pi_inv[v];};
            parlay::integer_sort_inplace(B[i],sort_f);
            T[i] = B[i][0];
        }
    });
    return parlay::reduce_max(bag_size);
}

}