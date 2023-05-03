#pragma once

namespace gbbs{

namespace pairing_heap{

template <typename key_type>
struct pairing_heap_node{
    key_type key;
    pairing_heap_node* child;
    pairing_heap_node* sibling;
    pairing_heap_node* parent;
    
    pairing_heap_node(key_type _key, pairing_heap_node* _child = nullptr, pairing_heap_node* _sibling = nullptr, 
        pairing_heap_node* _parent = nullptr) : key(_key), child(_child), sibling(_sibling), parent(_parent) {}
};

template <typename key_type>
pairing_heap_node<key_type>* meld(pairing_heap_node<key_type>* L, pairing_heap_node<key_type>* R){
    if (L == nullptr){
        return R;
    } else if (R == nullptr){
        return L;
    } else{
        if (L->key <= R->key){
            R->sibling = L->child;
            if (L->child){
                L->child->parent = R;
            }
            L->child = R;
            R->parent = L;
            return L;
        } else {
            L->sibling = R->child;
            if (R->child){
                R->child->parent = L;
            }
            R->child = L;
            L->parent = R;
            return R;
        }
    }
}

template <typename key_type>
pairing_heap_node<key_type>* two_pass_merge(pairing_heap_node<key_type>* heap){
    if (heap == nullptr){
        return heap;
    } else if (heap->sibling == nullptr){
        heap->parent = nullptr;
        return heap;
    } else{
        pairing_heap_node<key_type>* temp1 = heap->sibling;
        pairing_heap_node<key_type>* temp2 = heap->sibling->sibling;
        heap->sibling = nullptr;
        heap->parent = nullptr;
        temp1->sibling = nullptr;
        temp1->parent = nullptr;
        return meld(meld(heap, temp1), two_pass_merge(temp2));
    }
}

template <typename key_type>
struct pairing_heap{
    pairing_heap_node<key_type>* root;

    pairing_heap(pairing_heap_node<key_type>* _root = nullptr) : root(_root){}

    key_type find_min(){
        return root->key;
    }

    key_type delete_min(){
        if (root){
            pairing_heap_node<key_type>* temp = root;
            root = two_pass_merge(root->child);
            auto key = temp->key;
            delete temp;
            return key;
        } else{
            std::cout << "Error: Empty Heap" << std::endl;
            std::exit(1);
        }
    }

    bool is_empty(){
        return (root == nullptr);
    }

    pairing_heap_node<key_type>* emplace(key_type key){
        pairing_heap_node<key_type>* new_node = new pairing_heap_node<key_type>(key);
        root = meld(root, new_node);
        return new_node;
    }

    void merge(pairing_heap<key_type>* heap){
        root = meld(root, heap->root);
        heap->root = nullptr;
    }

    void decrease_key(pairing_heap_node<key_type>* node, key_type new_key){
        assert(node->key >= new_key);
        node->key = new_key;
        std::cout << "debug0" << std::endl;
        if (node->parent){
            if (node->sibling){
                node->sibling->parent = node->parent;
            }
            // std::cout << "debug1" << std::endl;
            if (node->parent->child == node){
                node->parent->child = node->sibling;
            } else {
                node->parent->sibling = node->sibling;
            }
            // std::cout << "debug2" << std::endl;
            node->parent = nullptr;
            node->sibling = nullptr;
            root = meld(root, node);
            // std::cout << "debug3" << std::endl;
        }
    }

    pairing_heap_node<key_type>* heapify(pairing_heap_node<key_type>** A, size_t n){
        if (n <= 64){
            auto temp = A[0];
            for (size_t i=0; i<n; i++){
                temp = meld(temp, A[i]);
            }
            return temp;
        } else {
            pairing_heap_node<key_type> *heap1, *heap2;
            parlay::par_do(
                [&](){heap1 = heapify(A, n/2);},
                [&](){heap2 = heapify(A + n/2, n - n/2);}
            );
            return meld(heap1, heap2);
        }
    }

    void create_heap(pairing_heap_node<key_type>** A, size_t n){
        root = heapify(A, n);
    }
};

} // namespace pairing_heap
} // namespace gbbs