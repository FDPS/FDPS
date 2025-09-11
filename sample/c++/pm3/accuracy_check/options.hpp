#pragma once
struct Params {
    //PS::S64 n_ptcl_glb = 16384;
    PS::S64 n_ptcl_glb = 4096;
    PS::F64 theta = 0.5;
    PS::S32 i_cut = 1;
    PS::S32 n_cell_1d = 8;
    //PS::F64 coef_ipg_size_limit = 0.5;
    PS::F64 negative_charge_ratio = 0.0;
    std::string input_file = "init.dat";
    bool read_file = false;
    void dump() {
        std::cout << "n_ptcl_glb= " << n_ptcl_glb << std::endl;
        std::cout << "theta= " << theta << std::endl;
        std::cout << "i_cut= " << i_cut << std::endl;
        std::cout << "n_cell_1d= " << n_cell_1d << std::endl;
        std::cout << "negative_charge_ratio= " << negative_charge_ratio << std::endl;
    }
};

void print_options() {
    std::cout << "  -h: show this help message and exit" << std::endl;
    std::cout << "  -n: total number of particles (default: 1024)" << std::endl;
    std::cout << "  -t: opening angle for the tree for PM3 (default: 0.5)" << std::endl;
    std::cout << "  -C: number of neighbor cells (default: 1)" << std::endl;
    std::cout << "  -c: total number of cells in the 1 dimensional direction (default: 8)" << std::endl;
    //std::cout << "  -i: input file name (option)" << std::endl;
    std::cout << "  -r: ratio of negative charged particles ((0.0, 1.0), default: 0.0)" << std::endl;
}

void select_options(Params &params, int argc, char **argv) {
    if (PS::Comm::getRank() == 0) {
        int opt;
        while ((opt = getopt(argc, argv, "hn:t:C:c:r:i:")) != -1) {
            switch (opt) {
                case 'h':
                    print_options();
                    exit(0);
                case 'n':
                    params.n_ptcl_glb = std::stol(optarg);
                    break;
                case 't':
                    params.theta = std::stod(optarg);
                    break;
                case 'C':
                    params.i_cut = std::stoi(optarg);
                    break;
                case 'c':
                    params.n_cell_1d = std::stoi(optarg);
                    break;
                case 'r':
                    params.negative_charge_ratio = std::stod(optarg);
                    break;
                case 'i':
                    params.input_file = std::string(optarg);
                    params.read_file = true;
                    break;
                default:
                    std::cerr << "Invalid option: " << opt << std::endl;
                    print_options();
                    exit(1);
            }
        }
    }
    PS::Comm::broadcast(&params, 1);
}
