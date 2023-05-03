#include "gbbs/gbbs.h"
#include "TreeDecomp/common/union.h"

namespace gbbs{

// MinDegree
template <class Graph>
sequence<uintE> mindegree(Graph& GA){
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  // if (G_copy == nullptr){
  auto G_copy = sequence<sequence<uintE>>(n);
  parallel_for(0, n, [&](uintE i){
    G_copy[i] = sequence<uintE>::uninitialized(GA.get_vertex(i).out_degree());
    auto map_f = [&](uintE u, uintE v, W w, size_t j){
        G_copy[i][j] = v;
    };
    GA.get_vertex(i).out_neighbors().map_with_index(map_f);
    parlay::integer_sort_inplace(G_copy[i]);
  });
  // }
  auto order = sequence<uintE>::uninitialized(n);
  auto degree = sequence<uintE>::uninitialized(n);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { degree[i] = GA.get_vertex(i).out_degree();});
  auto isnumbered = sequence<bool>::uninitialized(n);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { isnumbered[i] = 0;});

  for(size_t i=0; i<n; i++){
    auto next = parlay::min_element(degree)-degree.begin();
    order[i] = next;
    isnumbered[next] = 1;
    
    parallel_for(0, degree[next], [&](size_t j){
      auto v = G_copy[next][j];
      if (!isnumbered[v]){
        auto f = [&](uintE x){return 1;};
        auto l = [&](uintE x){
          if (x == next)
            return 0;
          return 1;
        };
        auto r = [&](uintE x){
          if (x == v)
            return 0;
          return 1;
        };
        G_copy[v] = union_flr_par(G_copy[v], G_copy[next], f, l, r);
        degree[v] = G_copy[v].size();
      }
    });
    degree[next] = UINT_E_MAX;
  }
  // auto bagsize = parlay::tabulate(n,[&](uintE i){return G_copy[i].size();});
  // auto mxtw = parlay::reduce_max(bagsize);
  return order;
}

}