#pragma once

#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "gbbs/gbbs.h"
#include "TreeDecomp/gavril.h"
#include "TreeDecomp/minfill.h"
#include "TreeDecomp/mcs.h"

namespace gbbs {

#define RMINDEG 0
#define RMINFILL 1
#define MCS 2
#define MCS_M 3
#define NDD 4
#define ALL 5

template <class Graph>
sequence<uintE> Ordering(Graph& GA, int order_heuristic = RMINDEG){
	if (order_heuristic == RMINDEG)
		return DegeneracyOrder(GA);
	else if (order_heuristic == RMINFILL)
		return minfill(GA);
	else if (order_heuristic == MCS)
		return mcs(GA);
	else if (order_heuristic == MCS_M)
		return mcs_m(GA);
	else if (order_heuristic == NDD)
		return DegeneracyOrder(GA);
	else
		return parlay::tabulate(GA.n,[&](uintE i){return i;});
}

void printOut(sequence<size_t> max_width, sequence<size_t> mean_width, sequence<size_t> median_width){
	std::cout << "### MINDEG: Max Bag Size = " << max_width[0] << ", Mean Bag Size = " << mean_width[0] 
	<< ", Median Bag Size = " << median_width[0] << std::endl;
	std::cout << "### MINFILL: Max Bag Size = " << max_width[1] << ", Mean Bag Size = " << mean_width[1] 
	<< ", Median Bag Size = " << median_width[1] << std::endl;
	std::cout << "### MCS: Max Bag Size = " << max_width[2] << ", Mean Bag Size = " << mean_width[2] 
	<< ", Median Bag Size = " << median_width[2] << std::endl;
	std::cout << "### MCS-M: Max Bag Size = " << max_width[3] << ", Mean Bag Size = " << mean_width[3] 
	<< ", Median Bag Size = " << median_width[3] << std::endl;
}

template <class Graph>
size_t TreeDecompEfficient(Graph& GA, int order_heuristic=ALL){
	using W = typename Graph::weight_type;
	size_t n = GA.n, treewidth;

	if (order_heuristic == ALL){
		auto max_width = sequence<size_t>(ALL);
		// auto mean_width = sequence<double>(ALL);
		// auto median_width = sequence<uintE>(ALL);
		parallel_for(0, ALL, [&](size_t heuristic){
			auto pi = Ordering(GA, heuristic);
			auto T = sequence<uintE>::uninitialized(n);
			auto B = sequence<sequence<uintE>>::uninitialized(n);
			size_t tw = GavrilTreeDecomp(GA, pi, B, T);
			max_width[heuristic] = tw;
			// mean_width[heuristic-1] = (double)parlay::reduce(bag_size)/n;
			// median_width[heuristic-1] = *parlay::kth_smallest(bag_size,n/2);
		});
		printOut(max_width, max_width, max_width); //TODO: median and mean width
		treewidth = parlay::reduce_min(max_width);
	}
	else{
		size_t max_width, mean_width, median_width;
		auto pi = Ordering(GA, order_heuristic);
		auto T = sequence<uintE>::uninitialized(n);
		auto B = sequence<sequence<uintE>>::uninitialized(n);
		size_t tw = GavrilTreeDecomp(GA, pi, B, T);
		max_width = tw;
		mean_width = tw; // (double)parlay::reduce(bag_size)/n;
		median_width = tw; // *parlay::kth_smallest(bag_size,n/2);

		std::cout << "### "<< order_heuristic << ": Max Bag Size = " << max_width << ", Mean Bag Size = " 
		<< mean_width << ", Median Bag Size = " << median_width << std::endl;
		treewidth = max_width;
	}
	return treewidth;
}

}  // namespace gbbs
