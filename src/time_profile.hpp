namespace ParticleSimulator {
struct TimeProfile {
    // for decompose domain
    F64 decompose_domain_all;
    F64 select_sample_particle;
    F64 decompose_domain;
    F64 decompose_domain_communication;
    F64 decompose_domain_sort_sample;

    // for exchange particle
    F64 exchange_particle;
    F64 exchange_particle_search;
    F64 exchange_particle_communication;

    // for interaction calculation
    F64 calc_force_all;
    F64 set_root_cell;
    F64 construct_lt;
    F64 set_particle_lt;
    F64 sort_particle_lt;
    F64 link_cell_lt;
    F64 set_moment_lt;
    F64 map_pm_cell_lt;
    F64 construct_exchange_let;
    F64 construct_let_1;
    F64 let_communication_1;
    F64 construct_let_2;
    F64 let_communication_2;
    F64 construct_gt;
    F64 set_particle_gt;
    F64 sort_particle_gt;
    F64 link_cell_gt;
    F64 set_moment_gt;
    F64 map_pm_cell_gt;
    F64 construct_ipg;
    F64 construct_interaction_list;
    F64 calculate_force;
    F64 write_back;

    //-----------------------
    // ParticleMeshMultipole
    //-----------------------
    F64 PMMM__calc_force_all;
    F64 PMMM__write_back;
    F64 PMMM__initialize;
    F64 PMMM__set_cell;
    F64 PMMM__set_ip_info_to_cell;
    F64 PMMM__calc_msum_and_quad0;
    F64 PMMM__calc_multipole_moment;
    F64 PMMM__collect_multipole_moment;
    F64 PMMM__GFC__initialize;
    F64 PMMM__GFC__calc_gf_r;
    F64 PMMM__GFC__redist_gf_r;
    F64 PMMM__GFC__calc_gf_k;
    F64 PMMM__GFC__redist_gf_k;
    F64 PMMM__M2L__initialize;
    F64 PMMM__M2L__preproc_mm_r_comm;   // redistMM
    F64 PMMM__M2L__redist_mm_r;         // redistMM
    F64 PMMM__M2L__postproc_mm_r_comm;  // redistMM
    F64 PMMM__M2L__mm_r_to_mm_k;
    F64 PMMM__M2L__gather_mm_k_trans;
    F64 PMMM__M2L__transform;
    F64 PMMM__M2L__scatter_le_k_trans;
    F64 PMMM__M2L__le_k_to_le_r;
    F64 PMMM__M2L__preproc_le_r_comm;   // redstLE
    F64 PMMM__M2L__redist_le_r;         // redistLE
    F64 PMMM__M2L__postproc_le_r_comm;  // redistLE
    F64 PMMM__dipole_correction;
    F64 PMMM__L2P;

    // F64 getTotalTime() const { return decompose_domain_all + exchange_particle + calc_force_all + write_back; }
    F64 getTotalTime() const { return decompose_domain_all + exchange_particle + calc_force_all + write_back; }

    void dump_domain_info(std::ostream &fout = std::cout, const S32 level = 0) const {
        fout << "decompose_domain_all: " << decompose_domain_all << std::endl;
        fout << "  select_sample_particle: " << select_sample_particle << std::endl;
        fout << "  decompose_domain: " << decompose_domain << std::endl;
        fout << "    communication: " << decompose_domain_communication << std::endl;
        fout << "    sort_sample: " << decompose_domain_sort_sample << std::endl;
    }

    void dump_particle_system(std::ostream &fout = std::cout, const S32 level = 0) const {
        fout << "exchange_particle: " << exchange_particle << std::endl;
        fout << "  exchange_particle_search: " << exchange_particle_search << std::endl;
        fout << "  exchange_particle_communication: " << exchange_particle_communication << std::endl;
    }

    void dump_tree_for_force(std::ostream &fout = std::cout, const S32 level = 0) const {
        fout << "calc_force_all: " << calc_force_all << std::endl;
        fout << "  set root cell: " << set_root_cell << std::endl;
        fout << "  LT construction: " << construct_lt << std::endl;
        fout << "    set particle LT: " << set_particle_lt << std::endl;
        fout << "    sort particle LT: " << sort_particle_lt << std::endl;
        fout << "    link cell LT: " << link_cell_lt << std::endl;
        fout << "    set moment LT: " << set_moment_lt << std::endl;
        fout << "    map PM cell LT: " << map_pm_cell_lt << std::endl;
        fout << "  LET construction exchange: " << construct_exchange_let << std::endl;
        fout << "    construct LET 1: " << construct_let_1 << std::endl;
        fout << "    LET communication 1: " << let_communication_1 << std::endl;
        fout << "    construct LET 2: " << construct_let_2 << std::endl;
        fout << "    LET communication 2: " << let_communication_2 << std::endl;
        fout << "  GT construction: " << construct_gt << std::endl;
        fout << "    set particle GT: " << set_particle_gt << std::endl;
        fout << "    sort particle GT: " << sort_particle_gt << std::endl;
        fout << "    link cell GT: " << link_cell_gt << std::endl;
        fout << "    set moment GT: " << set_moment_gt << std::endl;
        fout << "    map PM cell GT: " << map_pm_cell_gt << std::endl;
        fout << "  ip group construction: " << construct_ipg << std::endl;
        fout << "  interaction-list construction: " << construct_interaction_list << std::endl;
        fout << "  force calculation: " << calculate_force << std::endl;
        fout << "write_back: " << write_back << std::endl;
    }

    void dump(std::ostream &fout = std::cout, const S32 level = 0) const {
        fout << "total: " << getTotalTime() << std::endl;
        dump_domain_info(fout, level);
        dump_particle_system(fout, level);
        dump_tree_for_force(fout, level);
    }

    void clear() {
        decompose_domain_all = 0.0;
        select_sample_particle = 0.0;
        decompose_domain = 0.0;
        decompose_domain_communication = 0.0;
        decompose_domain_sort_sample = 0.0;

        exchange_particle = 0.0;
        exchange_particle_search = 0.0;
        exchange_particle_communication = 0.0;

        calc_force_all = 0.0;
        set_root_cell = 0.0;
        construct_lt = 0.0;
        set_particle_lt = 0.0;
        sort_particle_lt = 0.0;
        link_cell_lt = 0.0;
        set_moment_lt = 0.0;
        map_pm_cell_lt = 0.0;
        construct_exchange_let = 0.0;
        construct_let_1 = 0.0;
        let_communication_1 = 0.0;
        construct_let_2 = 0.0;
        let_communication_2 = 0.0;
        construct_gt = 0.0;
        set_particle_gt = 0.0;
        sort_particle_gt = 0.0;
        link_cell_gt = 0.0;
        set_moment_gt = 0.0;
        map_pm_cell_gt = 0.0;
        construct_ipg = 0.0;
        construct_interaction_list = 0.0;
        calculate_force = 0.0;
        write_back = 0.0;

        PMMM__calc_force_all = 0;
        PMMM__write_back = 0;
        PMMM__initialize = 0;
        PMMM__set_cell = 0;
        PMMM__set_ip_info_to_cell = 0;
        PMMM__calc_msum_and_quad0 = 0;
        PMMM__calc_multipole_moment = 0;
        PMMM__collect_multipole_moment = 0;
        PMMM__GFC__initialize = 0;
        PMMM__GFC__calc_gf_r = 0;
        PMMM__GFC__redist_gf_r = 0;
        PMMM__GFC__calc_gf_k = 0;
        PMMM__GFC__redist_gf_k = 0;
        PMMM__M2L__initialize = 0;
        PMMM__M2L__preproc_mm_r_comm = 0;
        PMMM__M2L__redist_mm_r = 0;
        PMMM__M2L__postproc_mm_r_comm = 0;
        PMMM__M2L__mm_r_to_mm_k = 0;
        PMMM__M2L__gather_mm_k_trans = 0;
        PMMM__M2L__transform = 0;
        PMMM__M2L__scatter_le_k_trans = 0;
        PMMM__M2L__le_k_to_le_r = 0;
        PMMM__M2L__preproc_le_r_comm = 0;
        PMMM__M2L__redist_le_r = 0;
        PMMM__M2L__postproc_le_r_comm = 0;
        PMMM__dipole_correction = 0;
        PMMM__L2P = 0;
    }

    TimeProfile() { clear(); }

    TimeProfile operator+(const TimeProfile &rhs) const {
        TimeProfile ret;
        ret.decompose_domain_all = this->decompose_domain_all + rhs.decompose_domain_all;
        ret.select_sample_particle = this->select_sample_particle + rhs.select_sample_particle;
        ret.decompose_domain = this->decompose_domain + rhs.decompose_domain;
        ret.decompose_domain_communication = this->decompose_domain_communication + rhs.decompose_domain_communication;
        ret.decompose_domain_sort_sample = this->decompose_domain_sort_sample + rhs.decompose_domain_sort_sample;

        ret.exchange_particle = this->exchange_particle + rhs.exchange_particle;
        ret.exchange_particle_search = this->exchange_particle_search + rhs.exchange_particle_search;
        ret.exchange_particle_communication = this->exchange_particle_communication + rhs.exchange_particle_communication;

        ret.calc_force_all = this->calc_force_all + rhs.calc_force_all;
        ret.set_root_cell = this->set_root_cell + rhs.set_root_cell;
        ret.construct_lt = this->construct_lt + rhs.construct_lt;
        ret.set_particle_lt = this->set_particle_lt + rhs.set_particle_lt;
        ret.sort_particle_lt = this->sort_particle_lt + rhs.sort_particle_lt;
        ret.link_cell_lt = this->link_cell_lt + rhs.link_cell_lt;
        ret.set_moment_lt = this->set_moment_lt + rhs.set_moment_lt;
        ret.map_pm_cell_lt = this->map_pm_cell_lt + rhs.map_pm_cell_lt;
        ret.construct_exchange_let = this->construct_exchange_let + rhs.construct_exchange_let;
        ret.construct_let_1 = this->construct_let_1 + rhs.construct_let_1;
        ret.let_communication_1 = this->let_communication_1 + rhs.let_communication_1;
        ret.construct_let_2 = this->construct_let_2 + rhs.construct_let_2;
        ret.let_communication_2 = this->let_communication_2 + rhs.let_communication_2;
        ret.construct_gt = this->construct_gt + rhs.construct_gt;
        ret.set_particle_gt = this->set_particle_gt + rhs.set_particle_gt;
        ret.sort_particle_gt = this->sort_particle_gt + rhs.sort_particle_gt;
        ret.link_cell_gt = this->link_cell_gt + rhs.link_cell_gt;
        ret.set_moment_gt = this->set_moment_gt + rhs.set_moment_gt;
        ret.map_pm_cell_gt = this->map_pm_cell_gt + rhs.map_pm_cell_gt;
        ret.construct_ipg = this->construct_ipg + rhs.construct_ipg;
        ret.construct_interaction_list = this->construct_interaction_list + rhs.construct_interaction_list;
        ret.calculate_force = this->calculate_force + rhs.calculate_force;
        ret.write_back = this->write_back + rhs.write_back;

        ret.PMMM__calc_force_all = this->PMMM__calc_force_all + rhs.PMMM__calc_force_all;
        ret.PMMM__write_back = this->PMMM__write_back + rhs.PMMM__write_back;
        ret.PMMM__initialize = this->PMMM__initialize + rhs.PMMM__initialize;
        ret.PMMM__set_cell = this->PMMM__set_cell + rhs.PMMM__set_cell;
        ret.PMMM__set_ip_info_to_cell = this->PMMM__set_ip_info_to_cell + rhs.PMMM__set_ip_info_to_cell;
        ret.PMMM__calc_msum_and_quad0 = this->PMMM__calc_msum_and_quad0 + rhs.PMMM__calc_msum_and_quad0;
        ret.PMMM__calc_multipole_moment = this->PMMM__calc_multipole_moment + rhs.PMMM__calc_multipole_moment;
        ret.PMMM__collect_multipole_moment = this->PMMM__collect_multipole_moment + rhs.PMMM__collect_multipole_moment;
        ret.PMMM__GFC__initialize = this->PMMM__GFC__initialize + rhs.PMMM__GFC__initialize;
        ret.PMMM__GFC__calc_gf_r = this->PMMM__GFC__calc_gf_r + rhs.PMMM__GFC__calc_gf_r;
        ret.PMMM__GFC__redist_gf_r = this->PMMM__GFC__redist_gf_r + rhs.PMMM__GFC__redist_gf_r;
        ret.PMMM__GFC__calc_gf_k = this->PMMM__GFC__calc_gf_k + rhs.PMMM__GFC__calc_gf_k;
        ret.PMMM__GFC__redist_gf_k = this->PMMM__GFC__redist_gf_k + rhs.PMMM__GFC__redist_gf_k;
        ret.PMMM__M2L__initialize = this->PMMM__M2L__initialize + rhs.PMMM__M2L__initialize;
        ret.PMMM__M2L__preproc_mm_r_comm = this->PMMM__M2L__preproc_mm_r_comm + rhs.PMMM__M2L__preproc_mm_r_comm;
        ret.PMMM__M2L__redist_mm_r = this->PMMM__M2L__redist_mm_r + rhs.PMMM__M2L__redist_mm_r;
        ret.PMMM__M2L__postproc_mm_r_comm = this->PMMM__M2L__postproc_mm_r_comm + rhs.PMMM__M2L__postproc_mm_r_comm;
        ret.PMMM__M2L__mm_r_to_mm_k = this->PMMM__M2L__mm_r_to_mm_k + rhs.PMMM__M2L__mm_r_to_mm_k;
        ret.PMMM__M2L__gather_mm_k_trans = this->PMMM__M2L__gather_mm_k_trans + rhs.PMMM__M2L__gather_mm_k_trans;
        ret.PMMM__M2L__transform = this->PMMM__M2L__transform + rhs.PMMM__M2L__transform;
        ret.PMMM__M2L__scatter_le_k_trans = this->PMMM__M2L__scatter_le_k_trans + rhs.PMMM__M2L__scatter_le_k_trans;
        ret.PMMM__M2L__le_k_to_le_r = this->PMMM__M2L__le_k_to_le_r + rhs.PMMM__M2L__le_k_to_le_r;
        ret.PMMM__M2L__preproc_le_r_comm = this->PMMM__M2L__preproc_le_r_comm + rhs.PMMM__M2L__preproc_le_r_comm;
        ret.PMMM__M2L__redist_le_r = this->PMMM__M2L__redist_le_r + rhs.PMMM__M2L__redist_le_r;
        ret.PMMM__M2L__postproc_le_r_comm = this->PMMM__M2L__postproc_le_r_comm + rhs.PMMM__M2L__postproc_le_r_comm;
        ret.PMMM__dipole_correction = this->PMMM__dipole_correction + rhs.PMMM__dipole_correction;
        ret.PMMM__L2P = this->PMMM__L2P + rhs.PMMM__L2P;

        return ret;
    }

    const TimeProfile &operator+=(const TimeProfile &rhs) {
        (*this) = (*this) + rhs;
        return (*this);
    }

    template <typename T>
    TimeProfile operator/(const T n) const {
        TimeProfile ret;
        ret.decompose_domain_all = this->decompose_domain_all / n;
        ret.select_sample_particle = this->select_sample_particle / n;
        ret.decompose_domain = this->decompose_domain / n;
        ret.decompose_domain_communication = this->decompose_domain_communication / n;
        ret.decompose_domain_sort_sample = this->decompose_domain_sort_sample / n;
        ret.exchange_particle = this->exchange_particle / n;
        ret.exchange_particle_search = this->exchange_particle_search / n;
        ret.exchange_particle_communication = this->exchange_particle_communication / n;
        ret.calc_force_all = this->calc_force_all / n;
        ret.set_root_cell = this->set_root_cell / n;
        ret.construct_lt = this->construct_lt / n;
        ret.set_particle_lt = this->set_particle_lt / n;
        ret.sort_particle_lt = this->sort_particle_lt / n;
        ret.link_cell_lt = this->link_cell_lt / n;
        ret.set_moment_lt = this->set_moment_lt / n;
        ret.map_pm_cell_lt = this->map_pm_cell_lt / n;
        ret.construct_exchange_let = this->construct_exchange_let / n;
        ret.construct_let_1 = this->construct_let_1 / n;
        ret.let_communication_1 = this->let_communication_1 / n;
        ret.construct_let_2 = this->construct_let_2 / n;
        ret.let_communication_2 = this->let_communication_2 / n;
        ret.construct_gt = this->construct_gt / n;
        ret.set_particle_gt = this->set_particle_gt / n;
        ret.sort_particle_gt = this->sort_particle_gt / n;
        ret.link_cell_gt = this->link_cell_gt / n;
        ret.set_moment_gt = this->set_moment_gt / n;
        ret.map_pm_cell_gt = this->map_pm_cell_gt / n;
        ret.construct_ipg = this->construct_ipg / n;
        ret.construct_interaction_list = this->construct_interaction_list / n;
        ret.calculate_force = this->calculate_force / n;
        ret.write_back = this->write_back / n;

        return ret;
    }

    template <typename T>
    TimeProfile operator*(const T n) const {
        TimeProfile ret;
        ret.decompose_domain_all = this->decompose_domain_all * n;
        ret.select_sample_particle = this->select_sample_particle * n;
        ret.decompose_domain = this->decompose_domain * n;
        ret.decompose_domain_communication = this->decompose_domain_communication * n;
        ret.decompose_domain_sort_sample = this->decompose_domain_sort_sample * n;
        ret.exchange_particle = this->exchange_particle * n;
        ret.exchange_particle_search = this->exchange_particle_search * n;
        ret.exchange_particle_communication = this->exchange_particle_communication * n;
        ret.calc_force_all = this->calc_force_all * n;
        ret.set_root_cell = this->set_root_cell * n;
        ret.construct_lt = this->construct_lt * n;
        ret.set_particle_lt = this->set_particle_lt * n;
        ret.sort_particle_lt = this->sort_particle_lt * n;
        ret.link_cell_lt = this->link_cell_lt * n;
        ret.set_moment_lt = this->set_moment_lt * n;
        ret.map_pm_cell_lt = this->map_pm_cell_lt * n;
        ret.construct_exchange_let = this->construct_exchange_let * n;
        ret.construct_let_1 = this->construct_let_1 * n;
        ret.let_communication_1 = this->let_communication_1 * n;
        ret.construct_let_2 = this->construct_let_2 * n;
        ret.let_communication_2 = this->let_communication_2 * n;
        ret.construct_gt = this->construct_gt * n;
        ret.set_particle_gt = this->set_particle_gt * n;
        ret.sort_particle_gt = this->sort_particle_gt * n;
        ret.link_cell_gt = this->link_cell_gt * n;
        ret.set_moment_gt = this->set_moment_gt * n;
        ret.map_pm_cell_gt = this->map_pm_cell_gt * n;
        ret.construct_ipg = this->construct_ipg * n;
        ret.construct_interaction_list = this->construct_interaction_list * n;
        ret.calculate_force = this->calculate_force * n;
        ret.write_back = this->write_back * n;

        return ret;
    }
};
}  // namespace ParticleSimulator