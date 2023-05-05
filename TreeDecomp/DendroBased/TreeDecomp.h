#pragma once

#include "gbbs/gbbs.h"
#include "TreeDecomp/common/utils.h"
#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "dendrogram.h"

namespace gbbs {

template <class Graph>
void TreeDecompDendro(Graph& GA, int order_heuristic=ALL){
	if (order_heuristic == ALL){
        for (int i=0; i<ALL; i++){
            TreeDecompDendro(GA, i);
        }
        return;
	}
    auto n = GA.n;
    timer t;
    t.start();
    // size_t max_width, median_width;
    // double mean_width;
    t.next("Init Time");

    // Step 1: Computing Elimination Order
    auto pi = Ordering(GA, order_heuristic);
    // for (size_t i=0; i<n; i++){
    //     std::cout << pi[i] << " ";
    // }
    // std::cout << std::endl;
    t.next("Ordering Time");

    // Step 2: Assign weights to edges
    reassign_weights(GA, pi);
    t.next("Reassigning Weights Time");

    // Step 3: Compute MST
    auto mst_edges = MinimumSpanningForest_boruvka::MinimumSpanningForest(GA);
    t.next("MST Time");

    // Step 4: Compute Dendrogram
    auto T = DendrogramSeqUF(mst_edges, n, pi);
    t.next("Dendro Time");
    for (size_t i=0; i< n; i++){
        std::cout << i << ": " << T[i] << std::endl;
    }
    // double tt = t.total_time();
    // max_width = tw;
    // mean_width = (double)parlay::reduce(bag_size_new)/n;
    // median_width = *parlay::kth_smallest(bag_size_new,n/2);
    // printOut(order_heuristic, max_width, mean_width, median_width, tt);
    return;
}

}  // namespace gbbs
