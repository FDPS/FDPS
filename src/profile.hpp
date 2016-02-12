#pragma once

#ifdef PROFILE_K
#include<fjcoll.h>
#endif
namespace PROFILE{
    template<class T>
    void Start(T x){
#ifdef PROFILE_K
        start_collection(x);
#endif
    }
    template<class T>
    void Stop(T x){
#ifdef PROFILE_K
        stop_collection(x);
#endif
    }

    class Profile{
        static const int SIZE = 1024;
    public:
        char prefix[SIZE];
        char collect_sample_particle[SIZE];
        char decompose_domain[SIZE];
        char exchange_particle[SIZE];
        char set_particle_local_tree[SIZE];
        char set_particle_global_tree[SIZE];
        char make_local_tree[SIZE];
        char make_global_tree[SIZE];
        char set_root_cell[SIZE];
        char calc_force[SIZE];
        char calc_moment_local_tree[SIZE];
        char calc_moment_global_tree[SIZE];
        char make_LET_1st[SIZE];
        char make_LET_2nd[SIZE];
        char exchange_LET_1st[SIZE];
        char exchange_LET_2nd[SIZE];
        char morton_sort_local_tree[SIZE];
        char link_cell_local_tree[SIZE];
        char morton_sort_global_tree[SIZE];
        char link_cell_global_tree[SIZE];
        void setPrefix(const char * str){
            sprintf(prefix, "%s", str);
#ifdef PROFILE_K
            sprintf(calc_force, "%s_calc_force", prefix);
#endif
        }

    };
}
