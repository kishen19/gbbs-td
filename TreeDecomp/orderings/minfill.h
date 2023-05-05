#include "TreeDecomp/common/union.h"

namespace gbbs{
namespace minfill{

template <class Graph, class Seq>
inline void CountTriangles(Graph& G, Seq* counts) {
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
	// std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
				// << " work per block = " << work_per_block << "\n";
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
inline size_t intersect_f_par(Seq* A, size_t nA, Seq* B, size_t nB, const F& f) {
  auto seqA = gbbs::make_slice<uintE>((uintE *)A, nA);
  auto seqB = gbbs::make_slice<uintE>((uintE *)B, nB);

  return intersection::merge(seqA, seqB, f);
}

template <class Graph>
sequence<uintE> minfill(Graph& GA){
  using W = typename Graph::weight_type;
  size_t n = GA.n;

  auto G_copy = sequence<sequence<uintE>>(n);
  parallel_for(0, n, [&](uintE i){
	G_copy[i] = sequence<uintE>::uninitialized(GA.get_vertex(i).out_degree());
	auto map_f = [&](uintE u, uintE v, W w, size_t j){
		G_copy[i][j] = v;
	};
	GA.get_vertex(i).out_neighbors().map_with_index(map_f);
	parlay::integer_sort_inplace(G_copy[i]); // Required for Union
  });

  // TODO: Lots of sequences, try to improve memory usage.
	auto order = sequence<uintE>::uninitialized(n); //Final Ordering
	auto degree = parlay::delayed_seq<intE>(n, [&](size_t i){return (intE)G_copy[i].size();}); //Degree value
	auto counts = sequence<intE>(n, 0); //triangles
	CountTriangles(GA, counts.begin()); 
	parallel_for(0, n, kDefaultGranularity, [&](size_t i){counts[i] = 2*counts[i]; }); // 2*num_triangles - easier to update during minfill
	auto isnumbered = sequence<bool>(n, 0);
	auto fill = parlay::delayed_seq<intE>(n, [&](size_t i){//Fill values
		if (isnumbered[i]){
			return INT_E_MAX;
		} else{
			return degree[i]*(degree[i]-1) - counts[i];
		}
	}); 
	auto false_f = [&](uintE x){return 0;};
	auto true_f = [&](uintE x){return 1;};
	// auto null_f = [&](uintE x){};
  for(size_t i=0; i<n; i++){
		auto next = parlay::min_element(fill) - fill.begin(); // Obtain minfill element
		order[i] = next; isnumbered[next] = 1;
		parallel_for(0, degree[next],[&](size_t j){
			auto v = G_copy[next][j];
			auto v_f = [&](uintE x){return (x==v? 0: 1);};
			auto new_nghs = union_flr_par(G_copy[v], G_copy[next], false_f, false_f, v_f);
			auto prev_nghs = union_flr_par(G_copy[v], G_copy[next], true_f, false_f, false_f);
			auto num_new_nghs = new_nghs.size();
			auto num_prev_nghs = prev_nghs.size();
			counts[v] = counts[v] - (2*num_prev_nghs); // remove triangles incident on next
			counts[v] = counts[v] + (num_new_nghs*(num_new_nghs-1) + 2*num_new_nghs*num_prev_nghs); // triangles with 2 new nghs + triangles with 1 new ngh
			parallel_for(0, degree[v], kDefaultGranularity, [&](size_t k){
				auto w = G_copy[v][k];
				if (w != next){
					if (parlay::find(G_copy[next], w) == G_copy[next].end()){
						auto temp = union_flr_par(new_nghs, G_copy[w], true_f, false_f, false_f).size();// TODO: Improve this step, just require num, not set
						gbbs::write_add(&counts[v], 2*temp);
						gbbs::write_add(&counts[w], temp);
					}
					else{
						auto temp = union_flr_par(prev_nghs, G_copy[w], true_f, false_f, false_f).size();// TODO: Improve this step, just require num, not set
						gbbs::write_add(&counts[v], num_prev_nghs - 1 - temp);
					}
				}
			});
		});
		parallel_for(0, degree[next], [&](size_t j){
			auto v = G_copy[next][j];
			auto next_f = [&](uintE x){ return (x==next ? 0 : 1);};
			auto v_f = [&](uintE x){return (x==v? 0: 1);};
			auto new_adjv = union_flr_par(G_copy[v], G_copy[next], true_f, next_f, v_f);
			G_copy[v] = new_adjv;
		});
  }
  return order;
}

}
}