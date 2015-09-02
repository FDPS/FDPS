#pragma once

#include<cstring>

namespace ParticleSimulator{
    
    class Timer{
    private:
        enum{
            SPLIT_CAPACITY = 128,
            STRING_SIZE = 1024,
        };
        double begin_time_;
        double end_time_;
        double split_time_[SPLIT_CAPACITY];
        double split_time_total_[SPLIT_CAPACITY];
        double total_time_;
        double total_time_sub_;
        char func_name_[SPLIT_CAPACITY][STRING_SIZE];
        size_t cnt_;
        size_t n_trial_;
    public:
        Timer(){
            begin_time_ = end_time_ = 0.0;
            for(S32 i=0; i<SPLIT_CAPACITY; i++){
                split_time_[i] = split_time_total_[i] = 0.0;
                func_name_[i][0] = '\0';
            }
            total_time_ = total_time_sub_ = 0.0;
            cnt_ = n_trial_ = 0; 
        }
        void initialize(std::ostream & fout){
            fout<<"Timer initialize"<<std::endl;
            begin_time_ = end_time_ = 0.0;
            for(S32 i=0; i<SPLIT_CAPACITY; i++){
                split_time_[i] = split_time_total_[i] = 0.0;
                func_name_[i][0] = '\0';
            }
            total_time_sub_ = 0.0;
            cnt_ = n_trial_ = 0; 
        }
        void reset(){ 
            for(S32 i=0; i<SPLIT_CAPACITY; i++){
                split_time_[i] = 0.0;
                func_name_[i][0] = '\0';
            }
            cnt_ = 0;
            n_trial_++;
        }
        void start(){ begin_time_ = GetWtime(); }
        void stop(){
            end_time_ = GetWtime();
            split_time_[cnt_++] = end_time_ - begin_time_;
        }
        void stop(const char * str){
            end_time_ = GetWtime();
            strcpy(func_name_[cnt_], str);
            split_time_[cnt_++] = end_time_ - begin_time_;
        }
        void restart(){
            stop();
            start();
        }
        void restart(const char * str){
            stop(str);
            start();
        }
        void dump(std::ostream & fout){
            double total_time_per_loop = 0.0;
            for(size_t i=0; i<cnt_; i++){
                total_time_per_loop += split_time_[i];
            }
            total_time_ += total_time_per_loop;
            total_time_sub_ += total_time_per_loop;
            double max_total_time_per_loop;
            int max_rank_total_time_per_loop;
            Comm::getMaxValue(total_time_per_loop, Comm::getRank(), max_total_time_per_loop, max_rank_total_time_per_loop);

            double max_total_time_sub;
            int max_rank_total_time_sub;
            Comm::getMaxValue(total_time_sub_, Comm::getRank(), max_total_time_sub, max_rank_total_time_sub);
            fout<<"total_time= "<<total_time_<<"   total_time_per_loop(ave)= "<<max_total_time_sub / n_trial_<<" (@"<<max_rank_total_time_per_loop
                <<")   max_total_time_per_loop= "<<max_total_time_per_loop<<" (@"<<max_rank_total_time_per_loop<<")"<<std::endl;

            double cumulative_time = 0.0;
            for(size_t i=0; i<cnt_; i++){

                double max_time = 0.0;
                int max_rank = 0;
                Comm::getMaxValue(split_time_[i], Comm::getRank(), max_time, max_rank);

                cumulative_time += split_time_[i];
                double max_time_cum = 0.0;
                int max_rank_cum = 0;
                Comm::getMaxValue(cumulative_time, Comm::getRank(), max_time_cum, max_rank_cum);

                split_time_total_[i] += split_time_[i];
                double max_time_split_total = 0.0;
                int max_rank_split_total = 0;
                Comm::getMaxValue(split_time_total_[i], Comm::getRank(), max_time_split_total, max_rank_split_total);
                if( func_name_[i][0] == '\0' ){
                    fout<<"split_time_["<<i<<"]= "<<split_time_[i]<<"   cum_time= "<<cumulative_time
                        <<"   fraction: "<<((double)split_time_[i])/total_time_per_loop<<"   max_time= "<<max_time<<" (@"<<max_rank<<")   max_time_cum= "
                        <<max_time_cum<<" (@"<<max_rank_cum<<")  "<<"   max_time_split_ave= "<<max_time_split_total/n_trial_<<" (@"<<max_rank_split_total<<")"<<std::endl;
                }
                else{
                    fout<<func_name_[i]<<": "<<split_time_[i]<<"   cum_time= "<<cumulative_time
                        <<"   fraction: "<<((double)split_time_[i])/total_time_per_loop<<"   max_time= "<<max_time<<" (@"<<max_rank<<")   max_time_cum= "
                        <<max_time_cum<<" (@"<<max_rank_cum<<")  "<<"   max_time_split_ave= "<<max_time_split_total/n_trial_<<" (@"<<max_rank_split_total<<")"<<std::endl;
                }
            }
        }
    };
}
