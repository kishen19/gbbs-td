#pragma once

#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "gbbs/gbbs.h"
#include "TreeDecomp/gavril.h"
#include "TreeDecomp/efficientTD.h"
#include "TreeDecomp/orderings/mindeg.h"
#include "TreeDecomp/orderings/minfill.h"
#include "TreeDecomp/orderings/mcs.h"

namespace gbbs {

#define DEGSORT 0
#define MINDEG 1
#define RMINDEG 2
#define MINFILL 3
#define RMINFILL 4
#define LEXBFS 5 // TODO
#define LEX_M 6 //TODO
#define MCS 7
#define MCS_M 8
#define NDD 9 // TODO
#define ALL 10

std::string names[ALL] = {"Sorted Degrees", "Mindeg", "R-Mindeg", "Minfill", "R-Minfill", "Lex-BFS", "Lex-M", "MCS", "MCS-M", "NDD"};

void printOut(sequence<size_t> max_width, sequence<double> mean_width, sequence<size_t> median_width, sequence<double> runtimes){
	std::cout << std::endl;
	for (auto i=0; i< ALL; i++){
		std::cout << "=> " << names[i] << ": Max Bag Size = " << max_width[i] << ", Mean Bag Size = " << mean_width[i] 
		<< ", Median Bag Size = " << median_width[i] << ", Time = " << runtimes[i] << std::endl;
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
	} else if (order_heuristic == MINDEG){
		return mindegree(GA);
	} else if (order_heuristic == RMINDEG){
		return DegeneracyOrder(GA);
	} else if (order_heuristic == MINFILL){ // TODO
		return minfill(GA);
	} else if (order_heuristic == RMINFILL){
		return rminfill(GA);
	} else if (order_heuristic == LEXBFS){ //TODO
		return mcs(GA);
	} else if (order_heuristic == LEX_M){ //TODO
		return mcs(GA);
	} else if (order_heuristic == MCS){
		return mcs(GA);
	} else if (order_heuristic == MCS_M){
		return mcs_m(GA);
	} else if (order_heuristic == NDD){ // TODO
		return DegeneracyOrder(GA);
	} else{
		return parlay::tabulate(n, [&](uintE i){return i;});
	}
}

template <class Graph>
size_t TreeDecompEfficient(Graph& GA, int order_heuristic=ALL){
	using W = typename Graph::weight_type;
	size_t n = GA.n, treewidth;

	if (order_heuristic == ALL){
		auto max_width = sequence<size_t>(ALL);
		auto mean_width = sequence<double>(ALL);
		auto median_width = sequence<size_t>(ALL);
		auto runtimes = sequence<double>(ALL);
		parallel_for(0, ALL, [&](size_t heuristic){
			auto T = sequence<uintE>::uninitialized(n);
			auto B = sequence<sequence<uintE>>(n);
			auto B_new = sequence<std::set<uintE>>(n);
			auto bag_size = sequence<size_t>::uninitialized(n);
			auto bag_size_new = sequence<size_t>::uninitialized(n);
			timer t;
			t.start();
			auto pi = Ordering(GA, heuristic);
			size_t tw = EfficientTreeDecomp(GA, pi, B, T, bag_size, B_new, bag_size_new);
			double tt = t.stop();
			max_width[heuristic] = tw;
			mean_width[heuristic] = (double)parlay::reduce(bag_size_new)/n;
			median_width[heuristic] = *parlay::kth_smallest(bag_size_new,n/2);
			runtimes[heuristic] = tt;
		});
		printOut(max_width, mean_width, median_width, runtimes);
		treewidth = parlay::reduce_min(max_width);
	}
	else{
		size_t max_width, median_width;
		double mean_width;
		auto T = sequence<uintE>::uninitialized(n);
		auto B = sequence<sequence<uintE>>(n);
		auto B_new = sequence<std::set<uintE>>(n);
		auto bag_size = sequence<size_t>::uninitialized(n);
		auto bag_size_new = sequence<size_t>::uninitialized(n);
		timer t;
		t.start();
		auto pi = Ordering(GA, order_heuristic);
		size_t tw = EfficientTreeDecomp(GA, pi, B, T, bag_size, B_new, bag_size_new);
		double tt = t.stop();
		max_width = tw;
		mean_width = (double)parlay::reduce(bag_size_new)/n;
		median_width = *parlay::kth_smallest(bag_size_new,n/2);

		std::cout << "=> "<< names[order_heuristic] << ": Max Bag Size = " << max_width << ", Mean Bag Size = " 
		<< mean_width << ", Median Bag Size = " << median_width << ", Time = " << tt << std::endl;
		treewidth = max_width;
	}
	return treewidth;
}

}  // namespace gbbs
