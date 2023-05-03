#pragma once

#include "gbbs/gbbs.h"
#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"

namespace gbbs {

template <class Seq>
sequence<uintE> DendrogramSeqUF(Seq& mst_edges, size_t& n, sequence<uintE>& pi){
	size_t m = mst_edges.size();
	auto pi_inv = sequence<uintE>::uninitialized(n);
    parallel_for(0, n, [&](size_t i) {pi_inv[pi[i]] = i;});
	
	timer t;
	t.start();

	// Step 1: Sorting the edges by (weight, index)
	auto indices = sequence<size_t>::from_function(m, [&](size_t i){return i;});
	auto comp = [&](const size_t& a, const size_t& b){
		auto w_a = std::get<2>(mst_edges[a]);
		auto w_b = std::get<2>(mst_edges[b]);
		if (w_a < w_b){
			return true;
		} else if (w_a == w_b){
			return (a < b);
		} else {
			return false;
		}
	};
	parlay::sort_inplace(indices, comp);
	t.next("Sorting Edges Time");

	// Step 2: Applying Union Find to the sorted sequence of edges
	auto uf = simple_union_find::SimpleUnionAsyncStruct(n);
	auto rep = sequence<uintE>::from_function(n, [&](uintE i){return i;});
	auto parents = sequence<uintE>::from_function(n, [&](size_t i){return i;}); //Output Parent Array
	auto heights = sequence<uintE>(n,0); // Heights of every node in the dendrogram

	for(size_t ind = 0; ind < m; ind++){
		auto i = indices[ind];
		uintE u,v;
		u = std::get<0>(mst_edges[i]);
		v = std::get<1>(mst_edges[i]);
		auto w = std::get<2>(mst_edges[i]);
		u = simple_union_find::find_compress(u, uf.parents);
		v = simple_union_find::find_compress(v, uf.parents);
		parents[rep[u]] = w; parents[rep[v]] = w;
		heights[u] = std::max(heights[u], heights[v])+1; heights[v] = heights[u];
		uf.unite(u,v);
		rep[u] = w; rep[v] = w;
	};
	t.next("UF Total Time");
	std::cout << std::endl << "=> Dendrogram Height = " << parlay::reduce_max(heights) << std::endl;
	return parents;
}

}  // namespace gbbs
