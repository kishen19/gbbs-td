#pragma once

#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "gbbs/gbbs.h"

namespace gbbs {

#define ALL 0
#define MINDEG 1
#define MINFILL 2
#define AMD 3
#define NDD 4

template <class Graph>
sequence<uintE> Ordering(Graph& GA, int order_heuristic = MINDEG){
	if (order_heuristic == MINDEG)
		return DegeneracyOrder(GA);
	else if (order_heuristic == MINFILL)
		return DegeneracyOrder(GA);
	else if (order_heuristic == AMD)
		return DegeneracyOrder(GA);
	else if (order_heuristic == NDD)
		return DegeneracyOrder(GA);
	else
		return parlay::tabulate(GA.n,[&](uintE i){return i;});
}

// Implements a parallel version of Gavril's algorithm for computing the
// treewidth of a graph.
template <class Graph>
uintE GavrilTreewidth(Graph& GA, int order_heuristic = ALL){
	using W = typename Graph::weight_type;
	size_t n = GA.n;
	uintE treewidth;

	if (order_heuristic == ALL){
		auto max_width = sequence<uintE>(4);
		auto mean_width = sequence<double>(4);
		auto median_width = sequence<uintE>(4);
		parallel_for(1, 5, [&](size_t heuristic){
			auto order = Ordering(GA, heuristic);
			auto vtx_to_position = sequence<uintE>(n);
			parallel_for(0, n, [&](size_t i) {
				uintE v = order[i];
				vtx_to_position[v] = i;
			});
			auto bag_size = sequence<size_t>(n);
			parallel_for(0, n, [&](size_t i) {
				uintE pos_u = vtx_to_position[i];
				auto vtx_f = [&](const uintE& u, const uintE& v, const W& wgh) {
					uintE pos_v = vtx_to_position[v];
					return pos_u < pos_v;
				};
				bag_size[pos_u] = GA.get_vertex(i).out_neighbors().count(vtx_f) + 1;
			});
			max_width[heuristic-1] = parlay::reduce_max(bag_size);
			mean_width[heuristic-1] = (double)parlay::reduce(bag_size)/n;
			median_width[heuristic-1] = *parlay::kth_smallest(bag_size,n/2);
		});
		std::cout << "### MINDEG: Max Bag Size = " << max_width[0] << ", Mean Bag Size = " << mean_width[0] 
		<< ", Median Bag Size = " << median_width[0] << std::endl;
		std::cout << "### MINFILL: Max Bag Size = " << max_width[1] << ", Mean Bag Size = " << mean_width[1] 
		<< ", Median Bag Size = " << median_width[1] << std::endl;
		std::cout << "### AMD: Max Bag Size = " << max_width[2] << ", Mean Bag Size = " << mean_width[2] 
		<< ", Median Bag Size = " << median_width[2] << std::endl;
		std::cout << "### NDD: Max Bag Size = " << max_width[3] << ", Mean Bag Size = " << mean_width[3] 
		<< ", Median Bag Size = " << median_width[3] << std::endl;
		treewidth = parlay::reduce_min(max_width);
		return treewidth;
	}
	else{
		auto order = Ordering(GA, order_heuristic);
		auto vtx_to_position = sequence<uintE>(n);
		parallel_for(0, n, [&](size_t i) {
			uintE v = order[i];
			vtx_to_position[v] = i;
		});
		auto bag_size = sequence<size_t>(n);
		parallel_for(0, n, [&](size_t i) {
			uintE pos_u = vtx_to_position[i];
			auto vtx_f = [&](const uintE& u, const uintE& v, const W& wgh) {
				uintE pos_v = vtx_to_position[v];
				return pos_u < pos_v;
			};
			bag_size[pos_u] = GA.get_vertex(i).out_neighbors().count(vtx_f) + 1;
		});
		auto max_width = parlay::reduce_max(bag_size);
		auto mean_width = (double)parlay::reduce(bag_size)/n;
		auto median_width = *parlay::kth_smallest(bag_size,n/2);
		std::cout << "### "<< order_heuristic << ": Max Bag Size = " << max_width << ", Mean Bag Size = " 
		<< mean_width << ", Median Bag Size = " << median_width << std::endl;
		treewidth = max_width;
	}
	return treewidth;
}



// Implements a parallel version of Gavril's algorithm for computing the tree 
// decomposition of a graph.
template <class Graph>
uintE GavrilTreeDecomp(Graph& GA, int order_heuristic = MINDEG) {
	using W = typename Graph::weight_type;
	size_t n = GA.n;
	auto order = Ordering(GA, order_heuristic);
	auto vtx_to_position = sequence<uintE>(n);

	parallel_for(0, n, [&](size_t i) {
		uintE v = order[i];
		vtx_to_position[v] = i;
	});

	auto bag_size = sequence<size_t>(n);

	parallel_for(0, n, [&](size_t i) {
		uintE pos_u = vtx_to_position[i];
		auto vtx_f = [&](const uintE& u, const uintE& v, const W& wgh) {
		uintE pos_v = vtx_to_position[v];
		return pos_u < pos_v;
		};
		bag_size[pos_u] = GA.get_vertex(i).out_neighbors().count(vtx_f) + 1;
	});
	// parallel_for(0, n, [&](size_t i) {
	//   uintE v = degeneracy_order[i];
	//   auto map_f = [&] (uintE u, uintE v, W w){
	//     vertices[(size_t)u].edges.push_back(edge{v,(C)0,(C)w,(C)0,nullptr});
	//   };
	//   G.get_vertex(u).out_neighbors().map(map_f,false);
	//   vtx_to_position[v] = i;
	// });

	uintE treewidth = parlay::reduce_max(bag_size);
	std::cout << "### Width of the Tree Decomposition is: " << treewidth-1 << std::endl;
	return treewidth-1;
}

}  // namespace gbbs
