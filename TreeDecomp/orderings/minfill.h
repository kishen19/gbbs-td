#include "gbbs/gbbs.h"

namespace gbbs{

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
			deg[v]--;
		};
		our_neighbors.map(map_f, true);  // run map sequentially
		counts[next] = n*n;
	}
	return order;
}

}