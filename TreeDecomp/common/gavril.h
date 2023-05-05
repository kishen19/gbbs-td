#pragma once

#include "TreeDecomp/common/td.h"

namespace gbbs{

template <class Graph>
size_t GavrilSimple(Graph& GA, sequence<uintE>& pi, sequence<sequence<uintE>>& B, sequence<uintE>& T, sequence<size_t>& bag_size){
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
            parlay::integer_sort_inplace(B[i], sort_f);
            T[i] = B[i][0];
        }
    });
    return parlay::reduce_max(bag_size);
}

template <class Edges>
TDEdges* GavrilEdges(Edges& edges, sequence<uintE>& pi, sequence<uintE>& pi_inv){
    using edge = std::pair<uintE, uintE>;
	size_t n = pi.size();
    size_t m = edges.size();

    auto sort_f = [&](const edge& e, const edge& f){
        if (e.first == f.first){
            return (pi_inv[e.second] < pi_inv[f.second]);
        } else{
            return (e.first < f.first);
        }
    };
    sort_inplace(edges, sort_f);

    auto TD = new TDEdges(n, m);
    TD->bag_size[edges[m-1].first] = m;
    parallel_for(0, m-1, [&](size_t i){
        assert(edges[i]!=edges[i+1]);
        if (edges[i].first != edges[i+1].first){
            TD->bag_size[edges[i].first] = i+1;
        }
    });
    TD->parent[edges[0].first] = edges[0].second;
    parallel_for(1, m, [&](size_t i){
        if (edges[i].first != edges[i-1].first){
            TD->bag_size[edges[i].first] -= i;
            TD->parent[edges[i].first] = edges[i].second;
        }
    });

    TD->edges = sequence<uintE>::from_function(m, [&](size_t i){return edges[i].second;});
    return TD;
}

} // namespace gbbs