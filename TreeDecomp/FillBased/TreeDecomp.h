#pragma once

#include "gbbs/gbbs.h"
#include "TreeDecomp/common/utils.h"
#include "TreeDecomp/common/gavril.h"
#include "fill_in.h"

namespace gbbs {

template <class Graph>
void TreeDecompFill(Graph& GA, int order_heuristic=ALL){
	if (order_heuristic == ALL){
        for (int i=0; i<ALL; i++){
            TreeDecompFill(GA, i);
        }
        return;
	}
    auto n = GA.n;
    timer t;
    t.start();
    size_t max_width, median_width;
    double mean_width;
    t.next("Init Time");

    // Step 1: Computing Elimination Order
    auto pi = Ordering(GA, order_heuristic);
    auto pi_inv = sequence<uintE>::uninitialized(n);
    parallel_for(0, n, [&](size_t i) {pi_inv[pi[i]] = i;});
    t.next("Ordering Time");

    // Step 2: Compute Fill Edges
    auto F_pi = fill_in(GA, pi, pi_inv);
    t.next("Fill-in Time");
    std::cout << "# Fill-in edges = " << F_pi.size() - GA.m << std::endl;

    // Step 3: Gavril on Triangulated Graph
    auto TD = GavrilEdges(F_pi, pi, pi_inv);
    t.next("Gavril Time");

    // double tt = t.total_time();
    // max_width = parlay::reduce_max(TD->bag_size);
    // mean_width = (double)parlay::reduce(TD->bag_size)/n;
    // median_width = *parlay::kth_smallest(TD->bag_size,n/2);
    for (size_t i=0; i< n; i++){
        std::cout << i << ": " << TD->parent[i] << std::endl;
    }
    // printOut(order_heuristic, max_width, mean_width, median_width, tt);
    return;
}

}  // namespace gbbs
