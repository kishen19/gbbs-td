#pragma once

#include "TreeDecomp/orderings/mindeg.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "TreeDecomp/orderings/minfill.h"
#include "TreeDecomp/orderings/lexbfs.h"
#include "TreeDecomp/orderings/mcs.h"

namespace gbbs{

#define DEGSORT 0
#define MINDEG 1
#define DEGENERACY 2
#define MINFILL 3
#define RMINFILL 4
#define LEXBFS 5
#define LEX_M 6
#define MCS 7
#define MCS_M 8
// #define NDD 9 // TODO
#define ALL 9

std::string names[ALL] = {"Sorted Deg", "Mindeg", "Degeneracy", "Minfill", "R-Minfill", "Lex-BFS", "Lex-M", "MCS", "MCS-M"};//, "NDD"};

inline void printOut(int& heuristic, size_t& max_width, double& mean_width, size_t& median_width, double& runtime){
	std::cout << "=> " << names[heuristic] << ": Max Bag Size = " << max_width << ", Mean Bag Size = " << mean_width
	<< ", Median Bag Size = " << median_width << ", Time = " << runtime << std::endl;
}

template <class Graph>
void reassign_weights(Graph& G, sequence<uintE>& pi_inv){
	auto n = G.n;
	auto f = [&](const uintE& u, const uintE& v){
		return std::make_tuple(v, std::max(pi_inv[u], pi_inv[v]));
	};
	parallel_for(0, n, [&](size_t i){
		G.get_vertex(i).out_neighbors().map_weights(f);
	});
}

template <class Graph>
sequence<uintE> Ordering(Graph& GA, int order_heuristic = DEGSORT){
	auto n = GA.n;

	if (order_heuristic == DEGSORT){
		auto sort_f = [&](uintE v){return GA.get_vertex(v).out_degree();};
		auto order = parlay::tabulate(n, [&](uintE i){return i;});
		parlay::integer_sort_inplace(order, sort_f);
		return order;
	} else if (order_heuristic == MINDEG){
		return mindeg::mindeg(GA);
	} else if (order_heuristic == DEGENERACY){
		return DegeneracyOrder(GA);
	} else if (order_heuristic == MINFILL){
		return minfill::minfill(GA);
	} else if (order_heuristic == RMINFILL){
		return minfill::rminfill(GA);
	} else if (order_heuristic == LEXBFS){
		return lex::lexbfs(GA);
	} else if (order_heuristic == LEX_M){
		return lex::lex_m(GA);
	} else if (order_heuristic == MCS){
		return mcs::mcs(GA);
	} else if (order_heuristic == MCS_M){
		return mcs::mcs_m(GA);
	// } else if (order_heuristic == NDD){ // TODO
	// 	return DegeneracyOrder(GA);
	} else{
		std::cout << "Invalid Heuristic" << std::endl;
		std::exit(1);
	}
}

}