#ifndef ITERATOR_LEVEL_H
#define ITERATOR_LEVEL_H

#include "iterator.h"


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

struct CacheBlock {
    int i_start, i_end;
    int j_start, j_end;
    int k_start, k_end;
    
    // Micro-wavefronts specifically for this block
    std::vector<std::vector<int>> micro_ijk_level; 
    
    // Flattened, block-local precomputed arrays to reduce memory streams
    // AoS (Array of Structures) to maximize cache line utilization
    struct NodeData {
        CUSTOMREAL fac_a, fac_b, fac_c, fac_f;
        CUSTOMREAL T0v, T0p, T0t, T0r;
        CUSTOMREAL fun, change;
    };
    std::vector<std::vector<NodeData>> micro_node_data;
    std::vector<std::vector<int>> micro_dump_ijk;
    std::vector<std::vector<int>> micro_dump_ip1, micro_dump_im1;
    std::vector<std::vector<int>> micro_dump_jp1, micro_dump_jm1;
    std::vector<std::vector<int>> micro_dump_kp1, micro_dump_km1;
};

class Iterator_level_1st_order_blocked : public Iterator_level {
public:
    Iterator_level_1st_order_blocked(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
    bool is_initialized = false;
    // THE MASTER SCHEDULE
    // Dimension 1: Sweep Direction (0 to 7)
    // Dimension 2: Macro-Level (Topological distance from the sweep origin)
    // Dimension 3: Independent CacheBlocks at this specific Macro-Level
    std::vector<std::vector<std::vector<CacheBlock>>> macro_levels_all_swp;

    // Initialization routine to build the schedule for all 8 directions
    void initialize_blocks(Grid& grid);

private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping

};

#endif // ITERATOR_LEVEL_H