#pragma once
#include<cassert>
namespace ParticleSimulator{
    template <typename T, int NMAX=200>
    class Stack{
    private:
        T val[NMAX];
        int num;
    public:
        //Stack() : num(0) {}
        void push(const T &v){
            assert(num < NMAX);
            val[num++] = v;
        }
        T pop(){
            return val[--num];
        }
        bool empty() const{
            return 0==num;
        }
        void init(){
            num = 0;
        }
    };
}
