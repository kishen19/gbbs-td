#pragma once

#include "gbbs/gbbs.h"
#include "TreeDecomp/common/gavril.h"
#include "TreeDecomp/common/utils.h"
#include "efficientTD.h"

namespace gbbs {

template <class Graph>
void TreeDecompRepair(Graph& GA, int order_heuristic=ALL){
	if (order_heuristic == ALL){
		for (int i=0; i<ALL; i++){
            TreeDecompRepair(GA, i);
        }
        return;
	}
	size_t n = GA.n;
	timer t;
	t.start();
	size_t max_width, median_width;
	double mean_width;
	auto T = sequence<uintE>::uninitialized(n);
	auto B = sequence<sequence<uintE>>(n);
	auto B_new = sequence<std::set<uintE>>(n);
	auto bag_size = sequence<size_t>::uninitialized(n);
	auto bag_size_new = sequence<size_t>::uninitialized(n);
	t.next("Init Time");
	auto pi = Ordering(GA, order_heuristic);
	t.next("Ordering Time");
	size_t tw = EfficientTreeDecomp(GA, pi, B, T, bag_size, B_new, bag_size_new);
	t.next("TD Time");
	double tt = t.total_time();
	max_width = tw;
	mean_width = (double)parlay::reduce(bag_size_new)/n;
	median_width = *parlay::kth_smallest(bag_size_new,n/2);
	// for (size_t i=0; i< n; i++){
    //     std::cout << i << ": " << T[i] << std::endl;
    // }
	// printOut(order_heuristic, max_width, mean_width, median_width, tt);
	return;
}

}  // namespace gbbs
