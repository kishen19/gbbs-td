#pragma once

namespace gbbs{
namespace mcs{

// Max Cardinality Search (MCS) Algorithm
struct queue_cell;

struct cell{
  cell* prev;
  cell* next;
  queue_cell* head;
  uintE v;
};

struct queue_cell{
  queue_cell* prev;
  queue_cell* next;
  cell* vertices;
  queue_cell(queue_cell* _prev, queue_cell* _next, cell* _vertices) 
    : prev(_prev), next(_next), vertices(_vertices){}
};

template <class Graph>
sequence<uintE> mcs(Graph& GA){
	using W = typename Graph::weight_type;
	size_t n = GA.n;

	auto order = sequence<uintE>::uninitialized(n);
  auto cells = sequence<cell>(n);
  queue_cell *queue = new queue_cell(nullptr, nullptr, &cells[0]);
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    if (i < n-1){
      cells[i].next = &cells[i+1];
    } else {
      cells[i].next = nullptr;
    }
    if (i > 0){
      cells[i].prev = &cells[i-1];
    } else {
      cells[i].prev = nullptr;
    }
    cells[i].head = queue;
    cells[i].v = i;
  });
  auto isnumbered = sequence<bool>(n,0);

  size_t ind = n;
	while (ind > 0){
    ind--;
		while (!queue->vertices){
      queue = queue->next;
      // delete queue->prev;
      queue->prev = nullptr;
    }
    auto next = queue->vertices->v;
    assert(isnumbered[next]==0);
    queue->vertices = (queue->vertices)->next;
    if (queue->vertices){
      (queue->vertices)->prev = nullptr;
    }
		order[ind] = next;
    isnumbered[next] = 1;
    
		auto map_f = [&](uintE u, uintE v, W wgh) {
      if (!isnumbered[v]){
        auto head = cells[v].head;
        if (!head->prev){ // if at top of queue, create a new queue cell
          queue_cell *new_queue_cell = new queue_cell(head->prev, head, nullptr);
          queue = new_queue_cell;
          head->prev = new_queue_cell;
        }
        // delete cells[v] from its set
        if (cells[v].prev){
          (cells[v].prev)->next = cells[v].next;
        } else{
          (cells[v].head)->vertices = cells[v].next;
        }
        if (cells[v].next){
          (cells[v].next)->prev = cells[v].prev;
        }
        // add cells[v] to new list
        cells[v].next = (head->prev)->vertices;
        if (cells[v].next){
          (cells[v].next)->prev = &cells[v];
        }
        cells[v].prev = nullptr;
        (head->prev)->vertices = &cells[v];
        cells[v].head = head->prev;
      }
		};
		GA.get_vertex(next).out_neighbors().map(map_f, false);
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
  // size_t reachable = 0;
  while (!Frontier.isEmpty()) {
    // reachable += Frontier.size();
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

// Max Cardinality Search (MCS) Algorithm - Old
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

} // namespace mcs
} // namespace gbbs