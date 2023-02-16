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
inline void CountTriangles(Graph& G, size_t* counts) {
  	using W = typename Graph::weight_type;
  	debug(std::cout << "Starting counting" << "\n";);
  	size_t n = G.n;

  	auto parallel_work = sequence<size_t>::uninitialized(n);
  	{
    	auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      		return G.get_vertex(v).out_degree();
    	};
		parallel_for(0, n, [&](size_t i) {
			auto monoid = parlay::plus<size_t>();
			parallel_work[i] = G.get_vertex(i).out_neighbors().reduce(map_f, monoid);
		});
  	}
	size_t total_work = parlay::scan_inplace(make_slice(parallel_work));

	size_t block_size = 50000;
	size_t n_blocks = total_work / block_size + 1;
	size_t work_per_block = (total_work + n_blocks - 1) / n_blocks;
	std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
				<< " work per block = " << work_per_block << "\n";
	auto f = [&](uintE u, uintE v, uintE w) {};
	auto run_intersection = [&](size_t start_ind, size_t end_ind) {
		for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
			auto our_neighbors = G.get_vertex(i).out_neighbors();
			size_t total_ct = 0;
			auto map_f = [&](uintE u, uintE v, W wgh) {
				auto their_neighbors = G.get_vertex(v).out_neighbors();
				total_ct += our_neighbors.intersect_f_par(&their_neighbors, f);
			};
			our_neighbors.map(map_f, false);  // run map sequentially
			counts[i] = total_ct/2;
		}
	};

	parallel_for(0, n_blocks, 1, [&](size_t i) {
		size_t start = i * work_per_block;
		size_t end = (i + 1) * work_per_block;
		auto less_fn = std::less<size_t>();
		size_t start_ind = parlay::binary_search(parallel_work, start, less_fn);
		size_t end_ind = parlay::binary_search(parallel_work, end, less_fn);
		run_intersection(start_ind, end_ind);
	});
}

template <class Graph>
sequence<uintE> minfill(Graph& GA){
	using W = typename Graph::weight_type;
	size_t n = GA.n;
	auto order = sequence<uintE>::uninitialized(n);
	auto deg = parlay::tabulate(n,[&](uintE i){
		return GA.get_vertex(i).out_degree();
	});
	auto counts = sequence<size_t>::uninitialized(n);
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) { counts[i] = 0; });
	CountTriangles(GA, counts.begin());
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) { counts[i] = deg[i]*(deg[i]-1)/2 - counts[i]; });
	size_t ind = 0;
	auto f = [&](uintE u, uintE v, uintE w) {};
	while (ind<n){
		auto next = parlay::min_element(counts)-counts.begin();
		order[ind] = next;
		ind++;
		auto our_neighbors = GA.get_vertex(next).out_neighbors();
		auto map_f = [&](uintE u, uintE v, W wgh) {
			auto their_neighbors = GA.get_vertex(v).out_neighbors();
			counts[v] += our_neighbors.intersect_f_par(&their_neighbors, f) + 1 - deg[v];
			deg[v] -= 1;
		};
		our_neighbors.map(map_f, true);  // run map sequentially
		counts[next] = n*n;
	}
	return order;
}

template <class Graph>
sequence<uintE> Ordering(Graph& GA, int order_heuristic = MINDEG){
	if (order_heuristic == MINDEG)
		return DegeneracyOrder(GA);
	else if (order_heuristic == MINFILL)
		return minfill(GA);
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
