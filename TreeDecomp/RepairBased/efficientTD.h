#pragma once

#include "gbbs/gbbs.h"
#include "TreeDecomp/common/gavril.h"
#include <set>

namespace gbbs{

template <class Graph>
size_t EfficientTreeDecomp(Graph& GA, sequence<uintE>& pi, sequence<sequence<uintE>>& B, sequence<uintE>& T, sequence<size_t>& bag_size, sequence<std::set<uintE>>& B_new, sequence<size_t>& bag_size_new){
    // using W = typename Graph::weight_type;
	size_t n = GA.n;
    timer t;
    t.start();
    GavrilSimple(GA, pi, B, T, bag_size);
    t.next("Gavril Time");
	auto pi_inv = sequence<uintE>::uninitialized(n+1);
    parallel_for(0, n, [&](size_t i) {pi_inv[pi[i]] = i;});
    pi_inv[n] = n;
    parallel_for(0, n, [&](size_t i) {bag_size_new[i] = bag_size[i];});
    parallel_for(0, n, [&](size_t i) {
        for (size_t j=0; j<bag_size[i]; j++){
            B_new[i].insert(B[i][j]);
        }
    });
    t.next("TD Init time");
    for (size_t i=0; i<n; i++){ //index i
        auto v = pi[i]; // vertex at index i
        if (bag_size_new[v]>0){
            for(auto it = B_new[v].begin(); it != B_new[v].end(); it++){
                auto cur = T[v];
                auto w = *it;
                while(cur != n && w!=n && w!= cur && B_new[cur].find(w)==B_new[cur].end()){
                    B_new[cur].insert(w);
                    bag_size_new[cur]++;
                    if (pi_inv[w] < pi_inv[T[cur]]){
                        auto temp = w;
                        w = T[cur];
                        T[cur] = temp;
                    }
                    cur = T[cur];
                }
            }
        }
    }
    t.next("TD Time");
    return parlay::reduce_max(bag_size_new);
}

}