#pragma once
namespace ParticleSimulator {

template <class Tforce, class Tepi, class Tepj, class Tmom = void, class Tsp = void>
class TreeForForceLong {
   public:
#if 1
    using Normal = TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, Tmom, Tmom, Tsp>;
    using WithCutoff = TreeForForce<SEARCH_MODE_LONG_CUTOFF, Tforce, Tepi, Tepj, Tmom, Tmom, Tsp>;
    using WithScatterSearch = TreeForForce<SEARCH_MODE_LONG_SCATTER, Tforce, Tepi, Tepj, Tmom, Tmom, Tsp>;
    using WithSymmetrySearch = TreeForForce<SEARCH_MODE_LONG_SYMMETRY, Tforce, Tepi, Tepj, Tmom, Tmom, Tsp>;
#else
    typedef TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, Tmom, Tmom, Tsp> Normal;
    typedef TreeForForce<SEARCH_MODE_LONG_CUTOFF, Tforce, Tepi, Tepj, Tmom, Tmom, Tsp> WithCutoff;
    typedef TreeForForce<SEARCH_MODE_LONG_SCATTER, Tforce, Tepi, Tepj, Tmom, Tmom, Tsp> WithScatterSearch;    // for P^3T
    typedef TreeForForce<SEARCH_MODE_LONG_SYMMETRY, Tforce, Tepi, Tepj, Tmom, Tmom, Tsp> WithSymmetrySearch;  // added by D.N. at 2018/11/28.
#endif
};

template <class Tforce, class Tepi, class Tepj>
class TreeForForceLong<Tforce, Tepi, Tepj, void, void> {
   public:
#if 1
    using MonopoleWithScatterSearch = TreeForForce<SEARCH_MODE_LONG_SCATTER, Tforce, Tepi, Tepj, MomentMonopole, MomentMonopole, SPJMonopole>;
    using QuadrupoleWithScatterSearch = TreeForForce<SEARCH_MODE_LONG_SCATTER, Tforce, Tepi, Tepj, MomentQuadrupole, MomentQuadrupole, SPJQuadrupole>;
    using MonopoleWithSymmetrySearch = TreeForForce<SEARCH_MODE_LONG_SYMMETRY, Tforce, Tepi, Tepj, MomentMonopole, MomentMonopole, SPJMonopole>;
    using QuadrupoleWithSymmetrySearch = TreeForForce<SEARCH_MODE_LONG_SYMMETRY, Tforce, Tepi, Tepj, MomentQuadrupole, MomentQuadrupole, SPJQuadrupole>;
    using Monopole = TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentMonopole, MomentMonopole, SPJMonopole>;
    using MonopoleWithCutoff = TreeForForce<SEARCH_MODE_LONG_CUTOFF, Tforce, Tepi, Tepj, MomentMonopole, MomentMonopole, SPJMonopole>;
    using Quadrupole = TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentQuadrupole, MomentQuadrupole, SPJQuadrupole>;
    using MonopoleGeometricCenter = TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentMonopoleGeometricCenter, MomentMonopoleGeometricCenter, SPJMonopoleGeometricCenter>;
    using DipoleGeometricCenter = TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentDipoleGeometricCenter, MomentDipoleGeometricCenter, SPJDipoleGeometricCenter>;
    using QuadrupoleGeometricCenter = TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentQuadrupoleGeometricCenter, MomentQuadrupoleGeometricCenter, SPJQuadrupoleGeometricCenter>;

#else
    // for P^3T
    typedef TreeForForce<SEARCH_MODE_LONG_SCATTER, Tforce, Tepi, Tepj, MomentMonopole, MomentMonopole, SPJMonopole> MonopoleWithScatterSearch;
    typedef TreeForForce<SEARCH_MODE_LONG_SCATTER, Tforce, Tepi, Tepj, MomentQuadrupole, MomentQuadrupole, SPJQuadrupole> QuadrupoleWithScatterSearch;
    typedef TreeForForce<SEARCH_MODE_LONG_SYMMETRY, Tforce, Tepi, Tepj, MomentMonopole, MomentMonopole, SPJMonopole> MonopoleWithSymmetrySearch;
    typedef TreeForForce<SEARCH_MODE_LONG_SYMMETRY, Tforce, Tepi, Tepj, MomentQuadrupole, MomentQuadrupole, SPJQuadrupole> QuadrupoleWithSymmetrySearch;

    // for P^3T + PM
    // typedef TreeForForce
    //<SEARCH_MODE_LONG_CUTOFF_SCATTER,
    // Tforce, Tepi, Tepj,
    // MomentMonopoleCutoffScatter,
    // MomentMonopoleCutoffScatter,
    // SPJMonopoleCutoffScatter> MonopoleWithCutoffScatterSearch;

    typedef TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentMonopole, MomentMonopole, SPJMonopole> Monopole;
    typedef TreeForForce<SEARCH_MODE_LONG_CUTOFF, Tforce, Tepi, Tepj, MomentMonopole, MomentMonopole, SPJMonopole> MonopoleWithCutoff;
    typedef TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentQuadrupole, MomentQuadrupole, SPJQuadrupole> Quadrupole;
    typedef TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentMonopoleGeometricCenter, MomentMonopoleGeometricCenter, SPJMonopoleGeometricCenter> MonopoleGeometricCenter;

    typedef TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentDipoleGeometricCenter, MomentDipoleGeometricCenter, SPJDipoleGeometricCenter> DipoleGeometricCenter;

    typedef TreeForForce<SEARCH_MODE_LONG, Tforce, Tepi, Tepj, MomentQuadrupoleGeometricCenter, MomentQuadrupoleGeometricCenter, SPJQuadrupoleGeometricCenter> QuadrupoleGeometricCenter;
#endif
};

template <class Tforce, class Tepi, class Tepj>
class TreeForForceShort {
   public:
    typedef TreeForForce<SEARCH_MODE_SYMMETRY, Tforce, Tepi, Tepj, MomentShort, MomentShort, SuperParticleBase> Symmetry;

    typedef TreeForForce<SEARCH_MODE_GATHER, Tforce, Tepi, Tepj, MomentShort, MomentShort, SuperParticleBase> Gather;

    typedef TreeForForce<SEARCH_MODE_SCATTER, Tforce, Tepi, Tepj, MomentShort, MomentShort, SuperParticleBase> Scatter;
};

}  // namespace ParticleSimulator