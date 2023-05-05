#pragma once

namespace gbbs{

template <class Graph>
sequence<std::pair<uintE, uintE>> fill_in(Graph& G, sequence<uintE>& pi, sequence<uintE>& pi_inv){
    using W = typename Graph::weight_type;
    auto n = G.n;

    auto f = sequence<uintE>::uninitialized(n);
    auto index = sequence<size_t>::uninitialized(n);
    auto F = sequence<std::pair<uintE, uintE>>();
    uintE w;
    for (size_t i = 0; i < n; i++){
        w = pi[i];
        f[w] = w;
        index[w] = i;
        auto map_f = [&](const uintE& u, const uintE& v, const W& wgh){
            if (pi_inv[v] < i){
                uintE x = v;
                while (index[x] < i){
                    index[x] = i;
                    F.emplace_back(x, u);
                    x = f[x];
                }
                if (f[x] == x){
                    f[x] = w;
                }
            }
        };
        G.get_vertex(w).out_neighbors().map(map_f, false);
    }
    return F;
}

} // namespace gbbs