#pragma once

namespace gbbs{

// struct vertex_bag{
//     uintE v;
//     uintE parent;
//     sequence<uintE> bag;

//     vertex_bag(uintE& _v, uintE& _parent = NULL): v(_v), parent(_parent){}
// }

struct TDEdges{
    size_t n;
    size_t m;
    sequence<uintE> parent; // parent bag id; bag id is self for root(s)
    sequence<uintE> edges;
    sequence<size_t> bag_size;

    TDEdges(size_t& _n, size_t& _m): n(_n), m(_m){
        parent = sequence<uintE>::from_function(n, [&](size_t i){return n;});
        edges = sequence<uintE>::uninitialized(m);
        bag_size = sequence<size_t>(n,0);
    }
};

} // namespace gbbs