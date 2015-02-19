#pragma once

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
        char func_name_[SPLIT_CAPACITY][STRING_SIZE];
        size_t cnt_;
    public:
        void reset(){ 
            for(S32 i=0; i<SPLIT_CAPACITY; i++){
                split_time_[i] = 0.0;
                func_name_[i][0] = '\0';
            }
            cnt_ = 0; 
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
            double total_time = 0.0;
            for(size_t i=0; i<cnt_; i++){
                total_time += split_time_[i];
            }
            double cumulative_time = 0.0;
            for(size_t i=0; i<cnt_; i++){
                cumulative_time += split_time_[i];
                double max_time = 0.0;
                int max_rank = 0;
                Comm::getMaxValue(split_time_[i], Comm::getRank(), max_time, max_rank);
                double max_time_cum = 0.0;
                int max_rank_cum = 0;
                Comm::getMaxValue(cumulative_time, Comm::getRank(), max_time_cum, max_rank_cum);
                if( func_name_[i][0] == '\0' ){
                    fout<<"split_time_["<<i<<"]= "<<split_time_[i]<<" cum_time= "<<cumulative_time
                        <<" fraction: "<<((double)split_time_[i])/total_time<<" max_time= "<<max_time<<" (@"<<max_rank<<") max_time_cum= "
                        <<max_time_cum<<" (@"<<max_rank_cum<<")"<<std::endl;
                }
                else{
                    fout<<func_name_[i]<<": "<<split_time_[i]<<" cum_time= "<<cumulative_time
                        <<" fraction: "<<((double)split_time_[i])/total_time<<" max_time= "<<max_time<<" (@"<<max_rank<<") max_time_cum= "
                        <<max_time_cum<<" (@"<<max_rank_cum<<")"<<std::endl;
                }
            }
        }


    };
}
