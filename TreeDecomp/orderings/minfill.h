#include "gbbs/gbbs.h"
#include "TreeDecomp/ligra_light.h"

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
sequence<uintE> rminfill(Graph& GA){
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
		our_neighbors.map(map_f);
		counts[next] = n*n;
	}
	return order;
}

template <class Seq, class F>
inline size_t intersect_f_par(Seq* A, uintE a, size_t nA, Seq* B, uintE b, size_t nB, const F& f) {
  auto seqA = gbbs::make_slice<uintE>((uintE *)A, nA);
  auto seqB = gbbs::make_slice<uintE>((uintE *)B, nB);

  auto merge_f = [&](uintE ngh) { f(a, b, ngh); };
  return intersection::merge(seqA, seqB, merge_f);
}

template <class Graph>
sequence<uintE> minfill(Graph& GA){
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  auto G_copy = sequence<sequence<uintE>>::uninitialized(n);
  parallel_for(0, n, [&](uintE i){
    G_copy[i] = sequence<uintE>::uninitialized(GA.get_vertex(i).out_degree());
    auto map_f = [&](uintE u, uintE v, W w, size_t j){
        G_copy[i][j] = v;
    };
    GA.get_vertex(i).out_neighbors().map_with_index(map_f);
    parlay::integer_sort_inplace(G_copy[i]);
  });
  auto order = sequence<uintE>::uninitialized(n); //Final Ordering
  auto degree = sequence<size_t>::uninitialized(n); //Degree value
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) { degree[i] = GA.get_vertex(i).out_degree();});
	auto counts = sequence<size_t>::uninitialized(n); //triangles
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) { counts[i] = 0; });
	CountTriangles(GA, counts.begin());
	auto fill = sequence<size_t>::uninitialized(n); //Fill values
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) { fill[i] = degree[i]*(degree[i]-1)/2 - counts[i]; });
  auto isnumbered = sequence<bool>::uninitialized(n);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { isnumbered[i] = 0;});

  for(size_t i=0; i<n; i++){
    auto next = parlay::min_element(fill) - fill.begin();
    order[i] = next;
    isnumbered[next] = 1;
    
    parallel_for(0, degree[next], [&](size_t j){
      auto v = G_copy[next][j];
      if (!isnumbered[v]){
        G_copy[v] = merge_out(G_copy[v], v, G_copy[next], next);
				auto added_nghs = G_copy[v].size() - degree[v];
				auto prev_nghs = degree[next] - added_nghs;
				counts[v] -= prev_nghs;
				degree[v] += added_nghs;
      }
    });

		// TODO: Use Vertex Subset and Edgemap to parallely recalculate the triangle counts for the vertices
		// in the first and second neighborhood of next.
		auto f = [&](uintE u, uintE v, uintE w) {};
		auto frontier_f = [&](uintE v){
			counts[v] = 0;
			for (size_t j=0; j < degree[v]; j++){
				auto w = G_copy[v][j];
				counts[v] += intersect_f_par(&G_copy[v], v, degree[v], &G_copy[w], w, degree[w], f);
			}
			counts[v]/=2;
			fill[v] = degree[v]*(degree[v]-1)/2 - counts[v];
		};
		auto visited = parlay::tabulate(n, [&](uintE v){return 0;});
		visited[next] = 1;

		auto edge_f = [&] (uintE u, uintE v) -> bool { 
			if (gbbs::atomic_compare_and_swap(&visited[v], 0, 1)){
				frontier_f(v);
				return 1;
			}
			else{
				return 0;
			}
		};
    auto cond_f = [&] (uintE v) { return !visited[v];};
    auto frontier_map = ligra::edge_map(G_copy, G_copy, edge_f, cond_f);

    // do the BFS
    auto frontier = ligra::vertex_subset<uintE>(next);
		// Update neighbors of next
		frontier = frontier_map(frontier);
		// Update neighbors of neighbors of next
		frontier = frontier_map(frontier);
    fill[next] = UINT_E_MAX;
  }
  return order;
}

}