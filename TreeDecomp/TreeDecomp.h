#pragma once

#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "gbbs/gbbs.h"
#include "TreeDecomp/gavril.h"
#include "TreeDecomp/minfill.h"
#include "TreeDecomp/mcs.h"

namespace gbbs {

#define DEGSORT 0
#define RMINDEG 1
#define RMINFILL 2
#define MCS 3
#define MCS_M 4
#define NDD 5
#define ALL 6

void printOut(sequence<size_t> max_width, sequence<double> mean_width, sequence<size_t> median_width){
	std::string names[ALL] = {"Sorted Degrees", "R-Mindeg", "R-Minfill", "MCS", "MCS-M", "NDD"};
	std::cout << std::endl;
	for (auto i=0; i< ALL; i++){
		std::cout << "=> " << names[i] << ": Max Bag Size = " << max_width[i] << ", Mean Bag Size = " << mean_width[i] 
		<< ", Median Bag Size = " << median_width[i] << std::endl;
	}
	std::cout << std::endl;
}

template <class Graph>
sequence<uintE> Ordering(Graph& GA, int order_heuristic = RMINDEG){
	auto n = GA.n;

	if (order_heuristic == DEGSORT){
		auto sort_f = [&](uintE v){return GA.get_vertex(v).out_degree();};
		auto order = parlay::tabulate(n, [&](uintE i){return i;});
		parlay::integer_sort_inplace(order, sort_f);
		return order;
	} else if (order_heuristic == RMINDEG)
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
		return parlay::tabulate(n, [&](uintE i){return i;});
}

template <class Graph>
size_t TreeDecompEfficient(Graph& GA, int order_heuristic=ALL){
	using W = typename Graph::weight_type;
	size_t n = GA.n, treewidth;

	if (order_heuristic == ALL){
		auto max_width = sequence<size_t>(ALL);
		auto mean_width = sequence<double>(ALL);
		auto median_width = sequence<size_t>(ALL);
		parallel_for(0, ALL, [&](size_t heuristic){
			auto pi = Ordering(GA, heuristic);
			auto T = sequence<uintE>::uninitialized(n);
			auto B = sequence<sequence<uintE>>::uninitialized(n);
			auto bag_size = sequence<size_t>::uninitialized(n);
			size_t tw = GavrilTreeDecomp(GA, pi, B, T, bag_size);
			max_width[heuristic] = tw;
			mean_width[heuristic] = (double)parlay::reduce(bag_size)/n;
			median_width[heuristic] = *parlay::kth_smallest(bag_size,n/2);
		});
		printOut(max_width, mean_width, median_width);
		treewidth = parlay::reduce_min(max_width);
	}
	else{
		size_t max_width, mean_width, median_width;
		auto pi = Ordering(GA, order_heuristic);
		auto T = sequence<uintE>::uninitialized(n);
		auto B = sequence<sequence<uintE>>::uninitialized(n);
		auto bag_size = sequence<size_t>::uninitialized(n);
		size_t tw = GavrilTreeDecomp(GA, pi, B, T, bag_size);
		max_width = tw;
		mean_width = (double)parlay::reduce(bag_size)/n;
		median_width = *parlay::kth_smallest(bag_size,n/2);

		std::cout << "=> "<< order_heuristic << ": Max Bag Size = " << max_width << ", Mean Bag Size = " 
		<< mean_width << ", Median Bag Size = " << median_width << std::endl;
		treewidth = max_width;
	}
	return treewidth;
}

}  // namespace gbbs
