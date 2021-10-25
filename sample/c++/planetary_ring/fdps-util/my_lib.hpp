#pragma once
namespace MY_LIB{
    namespace CONSTANT{
	constexpr long double G_cgs = 6.67259e-8L; // cm^3/(s^2 g)
	constexpr long double c_cgs = 2.99792458e10L; //  cm / s
	constexpr long double pi    = 3.14159265358979323846L;
	constexpr long double mass_saturn_cgs = 5.688e29L; // g
	constexpr long double au_cgs = 1.49597871475e13L; //cm
    };
    namespace LITERAL{
	// angle
	long double operator"" _deg2rad(long double deg){
	    return deg * CONSTANT::pi / 180.0;
	}
	long double operator"" _rad2deg(long double rad){
	    return rad * 180.0 / CONSTANT::pi;
	}
	long double operator"" _deg(long double deg){
	    return deg * CONSTANT::pi / 180.0;
	}
	
	// length
	long double operator"" _cm2m(long double len){
	    return len * 0.01;
	}
	long double operator"" _m2cm(long double len){
	    return len * 100.0;
	}
	long double operator"" _cm2km(long double len){
	    return len * 1e-5;
	}
	long double operator"" _km2cm(long double len){
	    return len * 1e5;
	}
	long double operator"" _cm2au(long double len){
	    return len / 1.49597871475e13;
	}
	long double operator"" _au2cm(long double len){
	    return len * 1.49597871475e13;
	}
	long double operator"" _cm2pc(long double len){
	    return len / 3.08567782e18;
	}
	long double operator"" _pc2cm(long double len){
	    return len * 3.08567782e18;
	}

	// mass
	long double operator"" _g2kg(long double mass){
	    return mass / 1000.0;
	}
	long double operator"" _kg2g(long double mass){
	    return mass * 1000.0;
	}
	long double operator"" _g2msun(long double mass){
	    return mass / 1.9884e33;
	}
	long double operator"" _msun2g(long double mass){
	    return mass * 1.9884e33;
	}
	long double operator"" _g2u(long double mass){
	    return mass / 1.66053906660e-30;
	}
	long double operator"" _u2g(long double mass){
	    return mass * 1.66053906660e-30;
	}
	
	// time
	long double operator"" _sec2yr(long double time){
	    return time / 3.15582e7;
	}
	long double operator"" _yr2sec(long double time){
	    return time * 3.15582e7;
	}
    }
    
    class UnitManager{
	double len_unit_;
	double mass_unit_;
	double time_unit_;
	double dens_unit_;
	double vel_unit_;
	double G_unit_;
	void calcLenUnit(){
	    len_unit_ = cbrt(G_unit_*time_unit_*time_unit_*mass_unit_);
	}
	void calcMassUnit(){
	    mass_unit_ = len_unit_*len_unit_*len_unit_ / (G_unit_ * time_unit_ * time_unit_);
	}
	void calcTimeUnit(){
	    time_unit_ = sqrt(len_unit_*len_unit_*len_unit_ / G_unit_ / mass_unit_);
	}
	void calcGUnit(){
	    G_unit_ = len_unit_*len_unit_*len_unit_ / (time_unit_*time_unit_ * mass_unit_);
	}
	void calcDensUnit(){
	    dens_unit_ = mass_unit_ / (len_unit_*len_unit_*len_unit_);
	}
	void calcVelUnit(){
	    vel_unit_ = len_unit_ / time_unit_;
	}
    public:
	void initializeByLMG(const double l, const double m, const double g){
	    len_unit_  = l;
	    mass_unit_ = m;
	    G_unit_    = g;
	    calcTimeUnit();
	    calcDensUnit();
	    calcVelUnit();
	}
	void initializeByLTG(const double l, const double t, const double g){
	    len_unit_  = l;
	    time_unit_ = t;
	    G_unit_    = g;
	    calcMassUnit();
	    calcDensUnit();
	    calcVelUnit();
	}
	void initializeByMTG(const double m, const double t, const double g){
	    mass_unit_  = m;
	    time_unit_  = t;
	    G_unit_     = g;
	    calcLenUnit();
	    calcDensUnit();
	    calcVelUnit();
	}
	double getLenUnit() const {
	    return len_unit_;
	}
	double getMassUnit() const {
	    return mass_unit_;
	}
	double getTimeUnit() const {
	    return time_unit_;
	}
	double getGUnit() const {
	    return G_unit_;
	}
	double getDensUnit() const {
	    return dens_unit_;
	}
	double getVelUnit() const {
	    return vel_unit_;
	}
	void dump(std::ostream & fout=std::cerr){
	    fout<<"len_unit_= "<<len_unit_
		<<" mass_unit_= "<<mass_unit_
		<<" time_unit_= "<<time_unit_
		<<" dens_unit_= "<<dens_unit_
		<<" vel_unit_= "<<vel_unit_
		<<" G_unit_= "<<G_unit_
		<<std::endl;
	}
	double getLen(const PS::F64 l) const {
	    return l * len_unit_;
	}
	double getMass(const PS::F64 m) const {
	    return m * mass_unit_;
	}
	double getTime(const PS::F64 t) const {
	    return t * time_unit_;
	}
	double getDens(const PS::F64 d) const {
	    return d * dens_unit_;
	}
	template<typename Tvec>
	Tvec getVel(const Tvec v) const {
	    return v * vel_unit_;
	}
    };

    void GetUniqueID(int id[], const int n_dst, const int min, const int max, int seed=0){
        const auto n_src = max - min + 1;
        assert(n_src >= n_dst);
        std::mt19937 mt(seed);
        std::uniform_int_distribution<int> dist(min, max);
        bool flag[n_src] = {false};
        for(int i=0; i<n_dst; i++){
            int j;
            do{
                j = dist(mt);
            }while(flag[j-min]==true);
            id[i] = j;
            flag[j-min] = true;
        }
    }
    
}
