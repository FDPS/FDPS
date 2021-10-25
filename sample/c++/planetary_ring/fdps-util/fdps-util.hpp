#pragma once
#include<tuple>
#include<particle_simulator.hpp>

namespace FDPS_UTIL{

    template<typename Tout, typename Tin>
    Tout LexicalCastImpl(Tin in){
	Tout out;
	std::stringstream ss;
	ss << in;
	ss >> out;
	return out;
    }
    
    template<typename Tin>
    struct MyLexical{
	Tin in_;
	MyLexical(Tin in) : in_(in) {}
	template<typename Tout>
	operator Tout() const {
	    return LexicalCastImpl<Tout>(in_);
	}
	template<typename Tout>
	operator Tout() {
	    return LexicalCastImpl<Tout>(in_);
	}
	
	template<typename T=Tin, typename std::enable_if<PS::my_is_char<T>::value>::type* = nullptr>
	operator bool() const{
	    if( strcmp(in_, "true") == 0 )
		return true;
	    else if( strcmp(in_, "false") == 0)
		return false;
	    else{
		std::cerr<<in_<<" should be \"true\" of \"flase\"\n";
		assert(0);
	    }
	}
	template<typename T=Tin, typename std::enable_if<PS::my_is_char<T>::value>::type* = nullptr>
	operator bool(){
	    if( strcmp(in_, "true") == 0 )
		return true;
	    else if( strcmp(in_, "false") == 0)
		return false;
	    else{
		std::cerr<<in_<<" should be \"true\" of \"flase\"\n";
		assert(0);
	    }
	}

	template<typename T=Tin, typename std::enable_if<PS::my_is_string<T>::value>::type* = nullptr>
	operator bool() const{
	    if( in_ == "true")
         	return true;
	    else if( in_ == "false")
		return false;
	    else{
		std::cerr<<in_<<" should be \"true\" of \"flase\"\n";
		assert(0);
	    }
	}
	template<typename T=Tin, typename std::enable_if<PS::my_is_string<T>::value>::type* = nullptr>
	operator bool(){
	    if( in_ == "true")
		return true;
	    else if( in_ == "false")
		return false;
	    else{
		std::cerr<<in_<<" should be \"true\" of \"flase\"\n";
		assert(0);
	    }
	}
    };

    template<typename Tin>
    MyLexical<Tin> LexicalCast(Tin in){
	return MyLexical<Tin>(in);
    }

    enum class CLO_TYPE{
	HAS_DEFAULT,
	HAS_NO_DEFAULT,
	HAS_NO_VALUE,
    };
    
    struct ComandLineOption{
	std::string long_name_;
	std::string description_;
	std::string val_str_;
	CLO_TYPE clo_type_;
	std::string short_name_;
	bool has_short_name_;
	bool need_;
	ComandLineOption(const std::string & ln, const std::string & des, const std::string & val, const CLO_TYPE & clo, const std::string sn, const bool has_sn, const bool nd)
	    : long_name_(ln), description_(des), val_str_(val), clo_type_(clo), short_name_(sn), has_short_name_(has_sn), need_(nd){}
	void setShortName(std::string && sn){
	    short_name_ = sn;
	    has_short_name_ = true;
	}
	void print(std::ostream & fout=std::cerr) const {
	    if( clo_type_ == CLO_TYPE::HAS_DEFAULT && has_short_name_ == true){
		printf("--%-15s-%-5s%-10s [default: %s]\n", long_name_.c_str(), short_name_.c_str(), description_.c_str(), val_str_.c_str());
	    }	    
	    else if( clo_type_ == CLO_TYPE::HAS_DEFAULT && has_short_name_ == false){
		printf("--%-21s%-10s [default: %s]\n", long_name_.c_str(), description_.c_str(), val_str_.c_str());
	    }
	    else{
		printf("--%-21s%-10s\n", long_name_.c_str(), description_.c_str());
	    }
	}
	bool needValue(){
	    return need_;
	}
	bool matchName(const std::string & name){
	    bool ret = false;
	    if(long_name_ == name || (has_short_name_ && short_name_ == name)){
		ret = true;
	    }
	    return ret;
	}
	
    };

    class ComandLineOptionManager{
	int argc_;
	char **argv_;
	std::vector<ComandLineOption> options_;
	void printHelpImpl(std::ostream & fout=std::cerr){
	    fout<<"*** comand line option ***"<<std::endl;
	    for(const auto & x : options_){
		x.print();
	    }
	    fout<<"**************************"<<std::endl;
	}
	void printHelp(std::ostream & fout=std::cerr){
	    for(auto i=0; i<argc_; i++){
		std::string str = "";
		auto is_option = chopOptionName(str, i);
		if(!is_option) continue;
		if(str == "help" || str == "h"){
		    if(PS::Comm::getRank()==0){
			printHelpImpl();
		    }
		    PS::Comm::barrier();
		    exit(1);
		}
	    }
	}
	
	bool chopOptionName(std::string & str, const int i){
	    bool ret = false;
	    if(argv_[i][0] == '-'){
		ret = true;
		str = std::string(argv_[i]).substr(1);
		if(argv_[i][1] == '-'){
		    str = std::string(argv_[i]).substr(2);
		    if(argv_[i][2] == '-' || argv_[i][2] == '\0'){
			ret = false;
		    }
		}
	    }
	    return ret;
	}
	void checkNeedToSetVaule(){
	    for(auto & opt : options_){
		if(opt.needValue()){
		    bool is_ok = false;
		    for(auto i=0; i<argc_; i++){
			std::string name = "";
			if(chopOptionName(name, i) && opt.matchName(name)){
			    is_ok = true;
			}
		    }
		    if(!is_ok){
			if(PS::Comm::getRank()==0){
			    std::cerr<<"need to specify the value for "<<opt.long_name_<<std::endl;
			}
			PS::Comm::barrier();
			exit(1);
		    }
		}
	    }
	}
	
	int getIndexFromName(std::string & name){
	    for(size_t i=0; i<options_.size(); i++){
		if(options_[i].long_name_ == name || (options_[i].has_short_name_ == true && options_[i].short_name_ == name) ){
		    return i;
		}
	    }
	    return -1;
	}
    public:
	ComandLineOptionManager(int ac, char *av[]) : argc_(ac), argv_(av) {}
	void appendImpl(const std::string & ln,
			const std::string & des,
			const std::string & val,
			const CLO_TYPE & clo,
			const std::string & sn,
			const bool has_sn,
			const bool nd){
	    options_.push_back(ComandLineOption(ln, des, val, clo, sn, has_sn, nd));
	}
	void appendFlag(const std::string & ln,
			const std::string & des){
	    appendImpl(ln, des, "false", CLO_TYPE::HAS_NO_VALUE, "", false, false);
	}
	void appendFlag(const std::string & ln,
			const std::string & sn,
			const std::string & des){
	    appendImpl(ln, des, "false", CLO_TYPE::HAS_NO_VALUE, sn, true, false);
	}
	void appendNoDefault(const std::string & ln,
			     const std::string & des,
			     const bool need){
	    appendImpl(ln, des, "", CLO_TYPE::HAS_NO_DEFAULT, "", false, need);
	}
	void appendNoDefault(const std::string & ln,
			     const std::string & sn,
			     const std::string & des,
			     const bool need){
	    appendImpl(ln, des, "", CLO_TYPE::HAS_NO_DEFAULT, sn, true, need);
	}
	void append(const std::string & ln,
		    const std::string & des,
		    const std::string & val){
	    appendImpl(ln, des, val, CLO_TYPE::HAS_DEFAULT, "", false, false);
	}
	void append(const std::string & ln,
		    const std::string & sn,
		    const std::string & des,
		    const std::string & val){
	    appendImpl(ln, des, val, CLO_TYPE::HAS_DEFAULT, sn, true, false);
	}
	
	void read(){
	    printHelp();
	    checkNeedToSetVaule();
	    for(auto i=0; i<argc_; i++){
		std::string str = "";
		auto is_option = chopOptionName(str, i);
		if(!is_option) continue;
		auto index = getIndexFromName(str);
		if(index >= 0){
		    if(options_[index].clo_type_ == CLO_TYPE::HAS_NO_VALUE){
			options_[index].val_str_ = "true";
		    } else if(i+1 < argc_){
			options_[index].val_str_ = std::string(argv_[i+1]);
		    }
		}
	    }
	}
	bool hasOption(const std::string & name){
	    bool ret = false;
	    for(auto i=0; i<argc_; i++){
		std::string str = "";
		auto is_option = chopOptionName(str, i);
		if(!is_option) continue;
		if(name == str){
		    ret = true;
		}
	    }
	    return ret;
	}

	MyLexical<std::string> get(const char * name){
	    std::string name_tmp = std::string(name);
	    auto index = getIndexFromName(name_tmp);
	    if(index >= 0){
		return MyLexical<std::string>(options_[index].val_str_);
	    }
	    else{
                std::cerr<<"name= "<<name<<std::endl;
		assert(0);
	    }
	}

	template<typename Tout>
	Tout get(const char * name){
	    std::string name_tmp = std::string(name);
	    auto index = getIndexFromName(name_tmp);
	    if(index >= 0){
		return LexicalCastImpl<Tout>(options_[index].val_str_);
	    }
	    else{
                std::cerr<<"name= "<<name<<std::endl;
		assert(0);
	    }
	}
	
    };


    template<typename... Args>
    class FileHeader{
	std::tuple<Args...> args_;
	bool flag_para;    
	template<int Nelm, int Ielm=0, typename std::enable_if<(Nelm > Ielm)>::type* = nullptr>
	void writeAsciiImpl(FILE* fp){
	    auto & x = std::get<Ielm>(args_);
	    x.writeAscii(fp);
	    writeAsciiImpl<Nelm, Ielm+1>(fp);
	}
	template<int Nelm, int Ielm=0, typename std::enable_if<(Nelm <= Ielm)>::type* = nullptr>
	void writeAsciiImpl(FILE* fp){}

        template<int Nelm, int Ielm=0, typename std::enable_if<(Nelm > Ielm)>::type* = nullptr>
	void writeBinaryImpl(FILE* fp){
	    auto & x = std::get<Ielm>(args_);
	    x.writeBinary(fp);
	    writeBinaryImpl<Nelm, Ielm+1>(fp);
	}
	template<int Nelm, int Ielm=0, typename std::enable_if<(Nelm <= Ielm)>::type* = nullptr>
	void writeBinaryImpl(FILE* fp){}

        template<int Nelm, int Ielm=0, typename std::enable_if<(Nelm > Ielm)>::type* = nullptr>
	void readAsciiImpl(FILE* fp){
	    auto & x = std::get<Ielm>(args_);
	    x.readAscii(fp);
	    readAsciiImpl<Nelm, Ielm+1>(fp);
	}
	template<int Nelm, int Ielm=0, typename std::enable_if<(Nelm <= Ielm)>::type* = nullptr>
	void readAsciiImpl(FILE* fp){}

        template<int Nelm, int Ielm=0, typename std::enable_if<(Nelm > Ielm)>::type* = nullptr>
	void readBinaryImpl(FILE* fp){
	    auto & x = std::get<Ielm>(args_);
	    x.readBinary(fp);
	    readBinaryImpl<Nelm, Ielm+1>(fp);
	}
	template<int Nelm, int Ielm=0, typename std::enable_if<(Nelm <= Ielm)>::type* = nullptr>
	void readBinaryImpl(FILE* fp){}	

public:
	FileHeader(Args & ... args) : args_(args...){}

	void writeAscii(FILE* fp){
	    constexpr auto Nelm = std::tuple_size<std::tuple<Args...>>::value;
	    writeAsciiImpl<Nelm>(fp);
	}
	void writeBinary(FILE* fp){
	    constexpr auto Nelm = std::tuple_size<std::tuple<Args...>>::value;
	    writeBinaryImpl<Nelm>(fp);
	}
	int readAscii(FILE* fp){
	    constexpr auto Nelm = std::tuple_size<std::tuple<Args...>>::value;
	    readAsciiImpl<Nelm>(fp);
	    if(flag_para)
		return std::get<0>(args_).getNumberOfParticleLocal();
	    else
		return std::get<0>(args_).getNumberOfParticleGlobal();
	}
	int readBinary(FILE* fp){
	    constexpr auto Nelm = std::tuple_size<std::tuple<Args...>>::value;
	    readBinaryImpl<Nelm>(fp);
	    if(flag_para)
		return std::get<0>(args_).getNumberOfParticleLocal();
	    else
		return std::get<0>(args_).getNumberOfParticleGlobal();
	}
	void setParallelIO(const bool f){
	    flag_para = f;
	}
	
	std::tuple<Args...> getArgs() const {
	    return args_;
	}
    };
    
    class SnapshotManager{
	std::string file_name_base; // from program
	PS::F64 dt_snp;   // from program
	PS::F64 time_snp; // from program or input file
	PS::S32 id_snp; // from program or input file
	bool flag_para;  // from program
	bool flag_bin;   // from program

	class HeaderParam{
	    PS::F64 time;
	    PS::S64 n_loc;
	    PS::S64 n_glb;
	    PS::S64 id_snp;
	public:
	    HeaderParam() : time(-1.0), n_loc(-1), n_glb(-1), id_snp(-1){}
	    HeaderParam(PS::F64 t, PS::S64 nl, PS::S64 ng, PS::S64 id) : time(t), n_loc(nl), n_glb(ng), id_snp(id){}
	    void writeAscii(FILE* fp){
		fprintf(fp, "%12.11e   %lld   %lld   %lld \n", time, n_loc, n_glb, id_snp);
	    }
	    void readAscii(FILE * fp){
		auto tmp = fscanf(fp, "%lf   %lld   %lld   %lld", &time, &n_loc, &n_glb, &id_snp);
		assert(tmp == 4);
	    }
	    void writeBinary(FILE * fp){
		fwrite(&time, sizeof(time), 1, fp);
		fwrite(&n_loc, sizeof(n_loc), 1, fp);
		fwrite(&n_glb, sizeof(n_glb), 1, fp);
		fwrite(&id_snp, sizeof(id_snp), 1, fp);
	    }
	    void readBinary(FILE * fp){
		size_t tmp = 0;
		tmp += fread(&time, sizeof(time), 1, fp);
		tmp += fread(&n_loc, sizeof(n_loc), 1, fp);
		tmp += fread(&n_glb, sizeof(n_glb), 1, fp);
		tmp += fread(&id_snp, sizeof(id_snp), 1, fp);
		assert(tmp == 4);
	    }
	    PS::S64 getIdSnp() const {
		return id_snp;
	    }
	    PS::S64 getNumberOfParticleLocal() const {
		return n_loc;
	    }
	    PS::S64 getNumberOfParticleGlobal() const {
		return n_glb;
	    }
	    PS::F64 getTimeSys() const {
		return time;
	    }
	};
	
    public:
	SnapshotManager(const std::string name, const PS::F64 dt, const bool flag_p=false, const bool flag_b=false) :
	    file_name_base(name), dt_snp(dt), flag_para(flag_p), flag_bin(flag_b){}
	void setIdSnp(const PS::S32 id){
	    id_snp = id;
	}
	void setTimeSnp(const PS::F64 time){
	    time_snp = time;
	}

	template<typename T>
	void setParams(const T & param){
	    id_snp   = param.getIdSnp();
	    time_snp = param.getTimeSnp();
	}
	
	template<typename Tsys, typename... Args>
	void write(Tsys & system, const PS::F64 time_sys, Args & ... args){
	    if(time_sys < time_snp) return;
	    std::ostringstream ss;
	    ss << std::setw(5) << std::setfill('0') << id_snp;
	    std::string file_name_prefix = file_name_base + ss.str();
	    HeaderParam header_param(time_sys, system.getNumberOfParticleLocal(), system.getNumberOfParticleGlobal(), id_snp);
	    FileHeader<HeaderParam, Args...> file_header(header_param, args...);
	    file_header.setParallelIO(flag_para);
	    if(flag_para){
		const char fmt[] = "%s_%05d_%05d.dat";
		if(flag_bin){
		    system.writeParticleBinary(file_name_prefix.c_str(), fmt, file_header);
		} else {
		    system.writeParticleAscii(file_name_prefix.c_str(), fmt, file_header);
		}
	    }
	    else{
		std::string file_name = file_name_prefix + ".dat";
		if(flag_bin){
		    system.writeParticleBinary(file_name.c_str(), file_header);
		} else {
		    system.writeParticleAscii(file_name.c_str(), file_header);
		}
	    }
	    time_snp += dt_snp;
	    id_snp++;
	}

	template<typename Tsys, typename... Args>
	void read(const std::string & read_file_name_base,
		  Tsys & system,
		  PS::F64 & time_sys,
		  Args & ... args){
	    HeaderParam header_param;
	    FileHeader<HeaderParam, Args...> file_header(header_param, args...);
	    file_header.setParallelIO(flag_para);
	    if(flag_para){
		const char fmt[] = "%s_%05d_%05d.dat";
		if(flag_bin){
		    system.readParticleBinary(read_file_name_base.c_str(), fmt, file_header);
		} else {
		    system.readParticleAscii(read_file_name_base.c_str(), fmt, file_header);
		}
	    }
	    else{
		if(flag_bin){
		    system.readParticleBinary(read_file_name_base.c_str(), file_header);
		} else {
		    system.readParticleAscii(read_file_name_base.c_str(), file_header);
		}
	    }
	    auto args_ret = std::tie(header_param, args...);
	    args_ret = file_header.getArgs();
	    id_snp   = header_param.getIdSnp() + 1;
	    time_sys = header_param.getTimeSys();
	    time_snp = time_sys + dt_snp;
	}
    };
}
