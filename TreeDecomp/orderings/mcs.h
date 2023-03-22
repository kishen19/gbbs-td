#include "gbbs/gbbs.h"
// #include "BFS/NonDeterministicBFS/BFS.h"

namespace gbbs{

// template <class W>
// struct BFS_F {
//   uintE* Parents;
//   uintE* weights;
//   BFS_F(uintE* _Parents) : Parents(_Parents) {}
// //   inline bool update(uintE s, uintE d, W w) {
// //     if (Parents[d] == UINT_E_MAX) {
// //       Parents[d] = s;
// //       return 1;
// //     } else {
// //       return 0;
// //     }
// //   }
//   inline bool updateAtomic(uintE s, uintE d, W w) {
// 	uintE newV, oldV;
// 	do {
// 		oldV = Parents[v];
// 		if (oldV != UINT_E_MAX && weights[s] >= weights[oldV])
// 			return 0;
// 	} while (!atomic_compare_and_swap(Parents[v], oldV, newV));
//     return 1;
//   }
//   inline bool cond(uintE d) { return (Parents[d] == UINT_E_MAX); }
// };

// template <class Graph>
// inline sequence<uintE> BFS(Graph& G, uintE src) {
//   using W = typename Graph::weight_type;
//   /* Creates Parents array, initialized to all -1, except for src. */
//   auto Parents =
//       sequence<uintE>::from_function(G.n, [&](size_t i) { return UINT_E_MAX; });
//   Parents[src] = src;

//   vertexSubset Frontier(G.n, src);
//   size_t reachable = 0;
//   while (!Frontier.isEmpty()) {
//     std::cout << Frontier.size() << "\n";
//     reachable += Frontier.size();
//     Frontier = edgeMap(G, Frontier, BFS_F<W>(Parents.begin()), -1,
//                        sparse_blocked | dense_parallel);
//   }
//   std::cout << "Reachable: " << reachable << "\n";
//   return Parents;
// }

template <class Graph>
sequence<uintE> mcs(Graph& GA){
	using W = typename Graph::weight_type;
	size_t n = GA.n;
	auto order = sequence<uintE>::uninitialized(n);
	auto w = sequence<intE>::uninitialized(n);
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) { w[i] = 0; });

    size_t ind = n;
	while (ind > 0){
        ind--;
		auto next = parlay::max_element(w)-w.begin();
		order[ind] = next;
		auto map_f = [&](uintE u, uintE v, W wgh) {
            if (w[v] >= 0)
			    w[v]++;
		};
		GA.get_vertex(next).out_neighbors().map(map_f, true);  // run map sequentially
		w[next] = -1;
	}
	return order;
}

template <class Graph>
sequence<uintE> mcs_m(Graph& GA){
	using W = typename Graph::weight_type;
	size_t n = GA.n;
	auto order = sequence<uintE>::uninitialized(n);
	auto w = sequence<intE>::uninitialized(n);
	parallel_for(0, n, kDefaultGranularity, [&](size_t i) { w[i] = 0; });

    size_t ind = n;
	while (ind > 0){
        ind--;
		auto next = parlay::max_element(w)-w.begin();
		order[ind] = next;
		auto map_f = [&](uintE u, uintE v, W wgh) {
            if (w[v] >= 0)
			    w[v]++;
		};
		GA.get_vertex(next).out_neighbors().map(map_f, true);  // run map sequentially
		w[next] = -1;
	}
	return order;
}

}