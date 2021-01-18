//
// Created by parthajit on 23/2/20.
//

#ifndef BPNET_BASIC_03_DISJOINTSET_H
class Disjoint_Set{
    int* _array;
    int _size;
public:
    Disjoint_Set(int size){
        this->_size = size;
        this->_array= new int[size];
        for (int i = 0; i < size; ++i) {
            _array[i] = i;
        }
    }
    void Union(int u, int v){
        Link(Find(u), Find(v));
    }
    int Find(int v){
        while(_array[v] != v){
            v = _array[v];
        }
        return v;
    }
    void Link(int u, int v){ // Links the two leaders.
        _array[v] = u;
    }
    ~Disjoint_Set(){
        delete [] _array;
        //cout<<" disjoint deleted";
    }
};
#define BPNET_BASIC_03_DISJOINTSET_H

#endif //BPNET_BASIC_03_DISJOINTSET_H
