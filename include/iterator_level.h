#ifndef ITERATOR_LEVEL_H
#define ITERATOR_LEVEL_H

#include "iterator.h"

#include "simd_conf.h"


class Iterator_level : public Iterator {
public:
    Iterator_level(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
protected:
    void do_sweep_adj(int, Grid&, InputParams&) override ; // do sweeping for adjoint routine
    virtual void do_sweep(int, Grid&, InputParams&) {}; // do sweeping
};

class Iterator_level_tele : public Iterator {
public:
    Iterator_level_tele(InputParams& , Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
protected:
    void do_sweep_adj(int, Grid&, InputParams&) override ; // do sweeping for adjoint routine
    virtual void do_sweep(int, Grid&, InputParams&) {}; // do sweeping
};

struct SIMDBlock {
    // Neighbor Indices for Gathering tau_loc
    int idx_c[NSIMD];
    int idx_ip1[NSIMD], idx_im1[NSIMD];
    int idx_jp1[NSIMD], idx_jm1[NSIMD];
    int idx_kp1[NSIMD], idx_km1[NSIMD];

    // Physical Parameters (Loaded contiguously in L1 Cache)
    CUSTOMREAL fac_a[NSIMD], fac_b[NSIMD], fac_c[NSIMD], fac_f[NSIMD];
    CUSTOMREAL T0v[NSIMD], T0p[NSIMD], T0t[NSIMD], T0r[NSIMD];
    CUSTOMREAL fun[NSIMD], change[NSIMD];
    
    // Tracks active lanes for boundary padding
    int valid_lanes; 
};

class Iterator_level_1st_order_opt : public Iterator_level {
public:
    // Dimension 1: 8 Sweep Directions
    // Dimension 2: Wavefront Levels
    // Dimension 3: SIMD Blocks (Iteration chunks)
    std::vector<std::vector<std::vector<SIMDBlock>>> vv_simd_blocks_all_swp;
    bool aosoa_built = false;
    Iterator_level_1st_order_opt(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
    void build_aosoa(Grid& grid);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
};

class Iterator_level_1st_order : public Iterator_level {
public:
    Iterator_level_1st_order(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
};

class Iterator_level_3rd_order : public Iterator_level {
public:
    Iterator_level_3rd_order(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
};

class Iterator_level_1st_order_upwind : public Iterator_level {
public:
    Iterator_level_1st_order_upwind(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
};

class Iterator_level_1st_order_tele : public Iterator_level_tele {
public:
    Iterator_level_1st_order_tele(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
};

class Iterator_level_3rd_order_tele : public Iterator_level_tele {
public:
    Iterator_level_3rd_order_tele(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
};

class Iterator_level_1st_order_upwind_tele : public Iterator_level_tele {
public:
    Iterator_level_1st_order_upwind_tele(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
};

#endif // ITERATOR_LEVEL_H