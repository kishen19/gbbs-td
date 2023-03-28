#include "gbbs/gbbs.h"

namespace gbbs{

constexpr const size_t _bs_merge_base = 32;
constexpr const size_t _seq_merge_thresh = 2048;


template<class Seq>
sequence<uintE> seq_merge_full(const Seq& A, const uintE& a, const Seq& B, const uintE& b) {
  size_t nA = A.size(), nB = B.size();
  auto C = sequence<uintE>::uninitialized(nA+nB);
  size_t i = 0, j = 0;
  size_t ct = 0;
  while (i < nA && j < nB) {
    const uintE& u = A[i];
    const uintE& v = B[j];
    if (u == b){
      i++;
    } else if (v == a){
      j++;
    } else if (u == v) {
      C[ct] = u;
      i++;
      j++;
      ct++;
    } else if (u < v) {
      C[ct] = u;
      i++; ct++;
    } else {
      C[ct] = v;
      j++; ct++;
    }
  }
  while (i < nA){
    if (A[i] != b){
      C[ct] = A[i];
      ct++; 
    }
    i++;
  }
  while (j < nB){
    if (B[j] != a){
      C[ct] = B[j];
      ct++;
    }
    j++;
  }
  auto output = parlay::tabulate(ct, [&](size_t i){return C[i];});
  return output;
}

template<class Seq>
sequence<uintE> seq_merge(const Seq& A, const uintE& a, const Seq& B, const uintE& b) {
  size_t nA = A.size(), nB = B.size();
  size_t ct = 0, j = 0;
  auto C = sequence<uintE>::uninitialized(nA+nB);
  bool flagA = 0, flagB = 0;
  for (size_t i = 0; i < nA; i++) {
    const uintE& u = A[i];
    if (u == b){
      flagA = 1;
      continue;
    }
    size_t mB = parlay::binary_search(B, u, std::less<uintE>());
    if (mB < nB){
      parallel_for(0, mB - j, [&](size_t k){
        if (!flagB && B[mB] > a){
          if (B[j+k]>a){
            C[ct+k-1] = B[j+k];
          } else if (B[j+k]<a){
            C[ct+k] = B[j+k];
          }
        } else{
          C[ct+k] = B[j+k];
        }
      });
      if (!flagB && B[mB] > a){
        flagB = 1;
        ct--;
      }
      ct += mB-j; j = mB;
      C[ct] = u;
      ct++;
      if (u == B[mB]){
        j++;
      }
    } else{
      parallel_for(0, mB - j, [&](size_t k){
        if (!flagB){
          if (B[j+k] > a){
            C[ct+k-1] = B[j+k];
          } else if (B[j+k] < a){
            C[ct+k] = B[j+k];
          }
        } else{
          C[ct+k] = B[j+k];
        }
      });
      ct += mB-j-(!flagB? 1 : 0); j = mB;
      parallel_for(0, nA - i, [&](size_t k){
        if (!flagA){
          if (A[i+k] > b){
            C[ct+k-1] = A[i+k];
          } else if (A[i+k] < b){
            C[ct+k] = A[i+k];
          }
        } else{
          C[ct+k] = A[i+k];
        }
      });
      ct += nA-i - (!flagA? 1: 0);
      break;
    }
  }
  parallel_for(0, nB-j, [&](size_t k){
    if (!flagB){
      if (B[j] > a){
        C[ct+k-1] = B[j+k];
      } else if (B[j] < a){
        C[ct+k] = B[j+k];
      }
    } else {
      C[ct+k] = B[j+k];
    }
  });
  auto output = parlay::tabulate(ct, [&](size_t i){return C[i];});
  return output;
}

template<class Seq>
sequence<uintE> merge_out(const Seq& A, uintE a, const Seq& B, uintE b){
  size_t nA = A.size();
  size_t nB = B.size();
  size_t nR = nA + nB;
  if (nR < _seq_merge_thresh) {  // handles (small, small) using linear-merge
    return seq_merge_full(A, a, B, b);
  } else if (nB < nA) {
    return merge_out(B, b, A, a);
  } else if (nA < _bs_merge_base){
    return seq_merge(A, a, B, b);
  } else {
    size_t mA = nA / 2;
    size_t mB = parlay::binary_search(B, A[mA], std::less<uintE>());
    sequence<uintE> m_left, m_right;
    par_do(
        [&]() {m_left = merge_out(A.cut(0, mA), a, B.cut(0, mB), b); },
        [&]() {m_right = merge_out(A.cut(mA, nA), a, B.cut(mB, nB), b);}
    );
    m_left.append(m_right);
    return m_left;
  }
}

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
        G_copy[v] = merge_out(G_copy[v], v, G_copy[next], next);
        degree[v] = G_copy[v].size();
      }
    });
    degree[next] = UINT_E_MAX;
  }
  auto bagsize = parlay::tabulate(n,[&](uintE i){return G_copy[i].size();});
  auto mxtw = parlay::reduce_max(bagsize);
  return order;
}

}