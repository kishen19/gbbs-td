#include "gbbs/gbbs.h"
#include "TreeDecomp/common/pairing_heap.h"

namespace gbbs{

// Max Cardinality Search (MCS) Algorithm
template <class Graph>
sequence<uintE> mcs_old(Graph& GA){
	using W = typename Graph::weight_type;
	size_t n = GA.n;
	auto order = sequence<uintE>::uninitialized(n);
	auto w = sequence<uintE>::uninitialized(n);
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) { w[i] = 1; });

  size_t ind = n;
	while (ind > 0){
    ind--;
		auto next = parlay::max_element(w)-w.begin();
		order[ind] = next;
		auto map_f = [&](uintE u, uintE v, W wgh) {
      if (w[v] > 0)
        w[v]++;
		};
		GA.get_vertex(next).out_neighbors().map(map_f);
		w[next] = 0;
	}
	return order;
}

// Max Cardinality Search (MCS) Algorithm
template <class Graph>
sequence<uintE> mcs(Graph& GA){
	using W = typename Graph::weight_type;
	size_t n = GA.n;

	auto order = sequence<uintE>::uninitialized(n);
	auto w = gbbs::new_array_no_init<pairing_heap::pairing_heap_node<uintE>*>(n);
  parallel_for(0, n, [&](size_t i){
    w[i] = new pairing_heap::pairing_heap_node<uintE>(n);
  });
  auto visited = sequence<bool>(n,0);
  auto heap = new pairing_heap::pairing_heap<uintE>();
  heap->create_heap(w, n);

  size_t ind = n;
	while (ind > 0){
    ind--;
		auto next = heap->delete_min();
		order[ind] = next;
		auto map_f = [&](uintE u, uintE v, W wgh) {
      if (visited[v] == 0)
        heap->decrease_key(w[v], w[v]->key-1);
		};
		GA.get_vertex(next).out_neighbors().map(map_f, false);
		visited[next] = 1;
	}
	return order;
}

// MCS-M Algorithm
template <class W>
struct MCS_M_BFS_F {
  uintE* weights; // max weight in the (best) path from each vertex to src
  uintE* mcs_weights; // Weights from the MCS procedure
  MCS_M_BFS_F(uintE* _weights, uintE* _mcs_weights) : weights(_weights), mcs_weights(_mcs_weights) {}
  inline bool update(uintE s, uintE d, W w){
    return updateAtomic(s, d, w);
  }
  inline bool updateAtomic(uintE s, uintE d, W w){
    uintE oldW, newW = std::max(weights[s], mcs_weights[s]);
    do {
      oldW = weights[d];
      if (oldW <= newW)
        return 0;
    } while (!gbbs::atomic_compare_and_swap(&weights[d], oldW, newW));
    return 1;
  }
  inline bool cond(uintE d) {return (mcs_weights[d] == UINT_E_MAX); }
};

template <class Graph>
inline void MCS_M_BFS(Graph& G, uintE src, sequence<uintE>& w) {
  using W = typename Graph::weight_type;
  /* Creates Weights array, initialized to all -1, except for src. */
  auto weights = sequence<uintE>::from_function(G.n, [&](size_t i) {return UINT_E_MAX; });
  weights[src] = 0;

  vertexSubset Frontier(G.n, src);
  size_t reachable = 0;
  while (!Frontier.isEmpty()) {
    reachable += Frontier.size();
    Frontier = edgeMap(G, Frontier, MCS_M_BFS_F<W>(weights.begin(), w.begin()), -1,
                       sparse_blocked | dense_parallel);
  }
  parallel_for(0, G.n, [&](size_t v){
    if (w[v] > weights[v])
      w[v]++;
  });
}

template <class Graph>
sequence<uintE> mcs_m(Graph& GA){
	// using W = typename Graph::weight_type;
	size_t n = GA.n;
	auto order = sequence<uintE>::uninitialized(n);
	auto w = sequence<uintE>::uninitialized(n);
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) { w[i] = 1; });
  size_t ind = n;
	while (ind > 0){
    ind--;
		auto next = parlay::max_element(w)-w.begin();
		order[ind] = next;
    MCS_M_BFS(GA, next, w);
		w[next] = 0;
	}
	return order;
}

}