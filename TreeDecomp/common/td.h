#pragma once

namespace gbbs{

struct vertex_bag{
    uintE v;
    uintE parent;
    sequence<uintE> bag;

    vertex_bag(uintE& _v, uintE& _parent = NULL): v(_v), parent(_parent){}
}

struct TreeDecomposition{
    size_t n;
    sequence<vertex_bag> bags;

    TreeDecomposition(size_t& _n): n(_n){
        bags = sequence<vertex_bag>::from_function(_n, [&](size_t i){
            auto 
        });
    }
}

}