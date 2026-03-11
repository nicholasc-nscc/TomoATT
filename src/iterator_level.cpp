#include "iterator_level.h"

#ifdef USE_SIMD
#include "vectorized_sweep.h"
#endif


Iterator_level::Iterator_level(InputParams& IP, Grid& grid, Source& src, IO_utils& io, const std::string& src_name, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                : Iterator(IP, grid, src, io, src_name, first_init, is_teleseismic_in, is_second_run_in) {
    // do nothing
}


void Iterator_level::do_sweep_adj(int iswp, Grid& grid, InputParams& IP){
    
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-1-kkr; 
            else            kkr = kkr;
            if (t_dirc < 0) jjt = nt-1-jjt; 
            else            jjt = jjt;  
            if (p_dirc < 0) iip = np-1-iip; 
            else            iip = iip;  

            if(iip < 0 || jjt < 0 || kkr < 0 || iip >= np || jjt >= nt || kkr >= nr) {
                std::cout << "ERROR: iip = " << iip << ", jjt = " << jjt << ", kkr = " << kkr << std::endl;
            }

            //
            // calculate stencils
            //
            if (iip != 0    && jjt != 0    && kkr != 0 \
                && iip != np-1 && jjt != nt-1 && kkr != nr-1) {
                // calculate stencils
                calculate_stencil_adj(grid, iip, jjt, kkr);
            } else {
                calculate_boundary_nodes_adj(grid, iip, jjt, kkr);
            }

        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

}


Iterator_level_tele::Iterator_level_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, const std::string& src_name, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                : Iterator(IP, grid, src, io, src_name, first_init, is_teleseismic_in, is_second_run_in) {
    // do nothing
}


void Iterator_level_tele::do_sweep_adj(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-1-kkr; 
            else            kkr = kkr;
            if (t_dirc < 0) jjt = nt-1-jjt; 
            else            jjt = jjt;  
            if (p_dirc < 0) iip = np-1-iip; 
            else            iip = iip;  

            //
            // calculate stencils
            //
            if (iip != 0    && jjt != 0    && kkr != 0 \
             && iip != np-1 && jjt != nt-1 && kkr != nr-1) {
                // calculate stencils
                calculate_stencil_adj(grid, iip, jjt, kkr);
            } else {
                calculate_boundary_nodes_adj(grid, iip, jjt, kkr);
            }
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

}

void Iterator_level_1st_order_blocked::initialize_blocks(Grid& grid) {
    const int BX = 16, BY = 16, BZ = 16;
    int nb_p = (np + BX - 1) / BX;
    int nb_t = (nt + BY - 1) / BY;
    int nb_r = (nr + BZ - 1) / BZ;

    macro_levels_all_swp.resize(8);

    for (int iswp = 0; iswp < 8; ++iswp) {
        int p_step = (iswp & 4) ? -1 : 1;
        int t_step = (iswp & 2) ? -1 : 1;
        int r_step = (iswp & 1) ? -1 : 1;
        int num_macro_levels = nb_p + nb_t + nb_r - 2;
        macro_levels_all_swp[iswp].resize(num_macro_levels);

        for (int I = 0; I < nb_p; ++I) {
            for (int J = 0; J < nb_t; ++J) {
                for (int K = 0; K < nb_r; ++K) {
                    int level_I = I;
                    int level_J = J;
                    int level_K = K;
                    
                    CacheBlock block;
                    block.i_start = I * BX; block.i_end = std::min((I + 1) * BX, np);
                    block.j_start = J * BY; block.j_end = std::min((J + 1) * BY, nt);
                    block.k_start = K * BZ; block.k_end = std::min((K + 1) * BZ, nr);
                    
                    int max_micro = (block.i_end - block.i_start) + (block.j_end - block.j_start) + (block.k_end - block.k_start) - 2;
                    
                    block.micro_valid_nodes.resize(max_micro, 0);
                    block.micro_dump_ijk.resize(max_micro); block.micro_dump_ip1.resize(max_micro); block.micro_dump_im1.resize(max_micro);
                    block.micro_dump_jp1.resize(max_micro); block.micro_dump_jm1.resize(max_micro);
                    block.micro_dump_kp1.resize(max_micro); block.micro_dump_km1.resize(max_micro);
                    
                    block.micro_fac_a.resize(max_micro); block.micro_fac_b.resize(max_micro); block.micro_fac_c.resize(max_micro); block.micro_fac_f.resize(max_micro);
                    block.micro_T0v.resize(max_micro); block.micro_T0p.resize(max_micro); block.micro_T0t.resize(max_micro); block.micro_T0r.resize(max_micro);
                    block.micro_fun.resize(max_micro); block.micro_change.resize(max_micro);
                    
                    macro_levels_all_swp[iswp][level_I + level_J + level_K].push_back(block);
                }
            }
        }
        
        int n_levels = ijk_for_this_subproc.size(); 
        for (int i_level = 0; i_level < n_levels; i_level++) {
            size_t n_nodes = ijk_for_this_subproc[i_level].size();            
            for (int i_node = 0; i_node < n_nodes; i_node++) {
                int iip, jjt, kkr;
                V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);
                int I = iip / BX; int J = jjt / BY; int K = kkr / BZ;
                int level_I = I;
                int level_J = J;
                int level_K = K;
                
                CacheBlock* t_blk = nullptr;
                for (auto& b : macro_levels_all_swp[iswp][level_I + level_J + level_K]) {
                    if (b.i_start == I * BX && b.j_start == J * BY && b.k_start == K * BZ) { t_blk = &b; break; }
                }
                
                if (t_blk) {
                    int m_i = iip - t_blk->i_start;
                    int m_j = jjt - t_blk->j_start;
                    int m_k = kkr - t_blk->k_start;
                    int m_lvl = m_i + m_j + m_k;
                    
                    t_blk->micro_dump_ijk[m_lvl].push_back(vv_i__j__k__.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_dump_ip1[m_lvl].push_back(vv_ip1j__k__.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_dump_im1[m_lvl].push_back(vv_im1j__k__.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_dump_jp1[m_lvl].push_back(vv_i__jp1k__.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_dump_jm1[m_lvl].push_back(vv_i__jm1k__.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_dump_kp1[m_lvl].push_back(vv_i__j__kp1.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_dump_km1[m_lvl].push_back(vv_i__j__km1.at(iswp).at(i_level)[i_node]);
                    
                    t_blk->micro_fac_a[m_lvl].push_back(vv_fac_a.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_fac_b[m_lvl].push_back(vv_fac_b.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_fac_c[m_lvl].push_back(vv_fac_c.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_fac_f[m_lvl].push_back(vv_fac_f.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_T0v[m_lvl].push_back(vv_T0v.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_T0p[m_lvl].push_back(vv_T0p.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_T0t[m_lvl].push_back(vv_T0t.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_T0r[m_lvl].push_back(vv_T0r.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_fun[m_lvl].push_back(vv_fun.at(iswp).at(i_level)[i_node]);
                    t_blk->micro_change[m_lvl].push_back(vv_change.at(iswp).at(i_level)[i_node]);
                }
            }
        }
        
        for (auto& macro_plane : macro_levels_all_swp[iswp]) {
            for (auto& b : macro_plane) {
                for (int m_lvl = 0; m_lvl < b.micro_dump_ijk.size(); m_lvl++) {
                    int curr_sz = b.micro_dump_ijk[m_lvl].size();
                    b.micro_valid_nodes[m_lvl] = curr_sz; 
                    if (curr_sz > 0 && curr_sz % NSIMD != 0) {
                        int pad_needed = NSIMD - (curr_sz % NSIMD);
                        for (int pad = 0; pad < pad_needed; pad++) {
                            b.micro_dump_ijk[m_lvl].push_back(b.micro_dump_ijk[m_lvl].back());
                            b.micro_dump_ip1[m_lvl].push_back(b.micro_dump_ip1[m_lvl].back());
                            b.micro_dump_im1[m_lvl].push_back(b.micro_dump_im1[m_lvl].back());
                            b.micro_dump_jp1[m_lvl].push_back(b.micro_dump_jp1[m_lvl].back());
                            b.micro_dump_jm1[m_lvl].push_back(b.micro_dump_jm1[m_lvl].back());
                            b.micro_dump_kp1[m_lvl].push_back(b.micro_dump_kp1[m_lvl].back());
                            b.micro_dump_km1[m_lvl].push_back(b.micro_dump_km1[m_lvl].back());
                            
                            b.micro_fac_a[m_lvl].push_back(b.micro_fac_a[m_lvl].back());
                            b.micro_fac_b[m_lvl].push_back(b.micro_fac_b[m_lvl].back());
                            b.micro_fac_c[m_lvl].push_back(b.micro_fac_c[m_lvl].back());
                            b.micro_fac_f[m_lvl].push_back(b.micro_fac_f[m_lvl].back());
                            b.micro_T0v[m_lvl].push_back(b.micro_T0v[m_lvl].back());
                            b.micro_T0p[m_lvl].push_back(b.micro_T0p[m_lvl].back());
                            b.micro_T0t[m_lvl].push_back(b.micro_T0t[m_lvl].back());
                            b.micro_T0r[m_lvl].push_back(b.micro_T0r[m_lvl].back());
                            b.micro_fun[m_lvl].push_back(b.micro_fun[m_lvl].back());
                            b.micro_change[m_lvl].push_back(b.micro_change[m_lvl].back());
                        }
                    }
                }
            }
        }
    }
}

// ERROR index!!!!!!!!!!!!!
Iterator_level_1st_order_blocked::Iterator_level_1st_order_blocked(InputParams& IP, Grid& grid, Source& src, IO_utils& io, const std::string& src_name, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level(IP, grid, src, io, src_name, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}

void Iterator_level_1st_order_blocked::do_sweep(int iswp, Grid& grid, InputParams& IP) {
    if (!use_gpu) {
        if (!is_initialized) { initialize_blocks(grid); is_initialized = true; }
        set_sweep_direction(iswp);

#if defined(USE_AVX512) || defined(USE_AVX)
        __mT v_DP_inv = _mmT_set1_pT(1.0/dp), v_DT_inv = _mmT_set1_pT(1.0/dt), v_DR_inv = _mmT_set1_pT(1.0/dr);
        __mT v_DP_inv_half = _mmT_set1_pT(1.0/dp*0.5), v_DT_inv_half = _mmT_set1_pT(1.0/dt*0.5), v_DR_inv_half = _mmT_set1_pT(1.0/dr*0.5);

        // MACRO: Safely buffers std::vector data into 64-byte aligned memory for _mmT_load_pT
        #define LOAD_SOA_ALIGNED(vec_name, array_name) \
            alignas(64) CUSTOMREAL local_##array_name[NSIMD]; \
            for(int l=0; l<NSIMD; l++) local_##array_name[l] = block.array_name[m_level][i_vec+l]; \
            __mT vec_name = _mmT_loadu_pT(local_##array_name);

        auto& macro_levels = macro_levels_all_swp[iswp];
        for (int b_level = 0; b_level < macro_levels.size(); b_level++) {
            
            #pragma omp parallel for schedule(dynamic)
            for (size_t b_idx = 0; b_idx < macro_levels[b_level].size(); b_idx++) {
                CacheBlock& block = macro_levels[b_level][b_idx];

                for (int m_level = 0; m_level < block.micro_dump_ijk.size(); m_level++) {
                    int valid_nodes = block.micro_valid_nodes[m_level];
                    if (valid_nodes == 0) continue;

                    int num_iter = block.micro_dump_ijk[m_level].size() / NSIMD;
                    for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {
                        int i_vec = _i_vec * NSIMD;

                        __mT v_c__ = load_mem_gen_to_mTd(grid.tau_loc, &block.micro_dump_ijk[m_level][i_vec]);
                        __mT v_p__ = load_mem_gen_to_mTd(grid.tau_loc, &block.micro_dump_ip1[m_level][i_vec]);
                        __mT v_m__ = load_mem_gen_to_mTd(grid.tau_loc, &block.micro_dump_im1[m_level][i_vec]);
                        __mT v__p_ = load_mem_gen_to_mTd(grid.tau_loc, &block.micro_dump_jp1[m_level][i_vec]);
                        __mT v__m_ = load_mem_gen_to_mTd(grid.tau_loc, &block.micro_dump_jm1[m_level][i_vec]);
                        __mT v___p = load_mem_gen_to_mTd(grid.tau_loc, &block.micro_dump_kp1[m_level][i_vec]);
                        __mT v___m = load_mem_gen_to_mTd(grid.tau_loc, &block.micro_dump_km1[m_level][i_vec]);

                        LOAD_SOA_ALIGNED(v_fac_a, micro_fac_a); LOAD_SOA_ALIGNED(v_fac_b, micro_fac_b); LOAD_SOA_ALIGNED(v_fac_c, micro_fac_c); LOAD_SOA_ALIGNED(v_fac_f, micro_fac_f);
                        LOAD_SOA_ALIGNED(v_T0v, micro_T0v); LOAD_SOA_ALIGNED(v_T0p, micro_T0p); LOAD_SOA_ALIGNED(v_T0t, micro_T0t); LOAD_SOA_ALIGNED(v_T0r, micro_T0r);
                        LOAD_SOA_ALIGNED(v_fun, micro_fun); LOAD_SOA_ALIGNED(v_change, micro_change);

                        __mT v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2;
                        vect_stencil_1st_pre_simd(v_c__, v_p__, v_m__, v__p_, v__m_, v___p, v___m,
                                                  v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2,
                                                  v_DP_inv, v_DT_inv, v_DR_inv, v_DP_inv_half, v_DT_inv_half, v_DR_inv_half,
                                                  loc_I, loc_J, loc_K);

                        vect_stencil_1st_3rd_apre_simd(v_c__, v_fac_a, v_fac_b, v_fac_c, v_fac_f,
                                                       v_T0v, v_T0p, v_T0t, v_T0r, v_fun, v_change,
                                                       v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2,
                                                       v_DP_inv, v_DT_inv, v_DR_inv);

                        alignas(64) CUSTOMREAL local_dump_c[NSIMD];
                        _mmT_store_pT(local_dump_c, v_c__);
                        for (int i = 0; i < NSIMD; i++) {
                            if(i_vec+i >= valid_nodes) break; 
                            grid.tau_loc[block.micro_dump_ijk[m_level][i_vec+i]] = local_dump_c[i];
                        }
                    }
                }
            }
            synchronize_all_sub(); 
        }
#endif
    } else {
        // GPU fallback
    }
    if (subdom_main) calculate_boundary_nodes(grid);
}


// ERROR index!!!!!!!!!!!!!
Iterator_level_1st_order::Iterator_level_1st_order(InputParams& IP, Grid& grid, Source& src, IO_utils& io, const std::string& src_name, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level(IP, grid, src, io, src_name, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_1st_order::do_sweep(int iswp, Grid& grid, InputParams& IP){

    if(!use_gpu){

#if !defined USE_SIMD

        // set sweep direction
        set_sweep_direction(iswp);

        int iip, jjt, kkr;
        int n_levels = ijk_for_this_subproc.size();

        for (int i_level = 0; i_level < n_levels; i_level++) {
            size_t n_nodes = ijk_for_this_subproc[i_level].size();

            for (size_t i_node = 0; i_node < n_nodes; i_node++) {

                V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

                if (r_dirc < 0) kkr = nr-kkr; //kk-1;
                else            kkr = kkr-1;  //nr-kk;
                if (t_dirc < 0) jjt = nt-jjt; //jj-1;
                else            jjt = jjt-1;  //nt-jj;
                if (p_dirc < 0) iip = np-iip; //ii-1;
                else            iip = iip-1;  //np-ii;

                //
                // calculate stencils
                //
                if (grid.is_changed[I2V(iip, jjt, kkr)]) {
                    calculate_stencil_1st_order(grid, iip, jjt, kkr);
                } // is_changed == true
            } // end ijk

            // mpi synchronization
            synchronize_all_sub();

        } // end loop i_level

#elif USE_AVX512 || USE_AVX

        // preload constants
        __mT v_DP_inv      = _mmT_set1_pT(1.0/dp);
        __mT v_DT_inv      = _mmT_set1_pT(1.0/dt);
        __mT v_DR_inv      = _mmT_set1_pT(1.0/dr);
        __mT v_DP_inv_half = _mmT_set1_pT(1.0/dp*0.5);
        __mT v_DT_inv_half = _mmT_set1_pT(1.0/dt*0.5);
        __mT v_DR_inv_half = _mmT_set1_pT(1.0/dr*0.5);

        // store stencil coefs
        __mT v_pp1;
        __mT v_pp2;
        __mT v_pt1;
        __mT v_pt2;
        __mT v_pr1;
        __mT v_pr2;

        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc.at(i_level).size();

            int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

            // make alias to preloaded data
            __mT* v_fac_a  = (__mT*) vv_fac_a.at(iswp).at(i_level);
            __mT* v_fac_b  = (__mT*) vv_fac_b.at(iswp).at(i_level);
            __mT* v_fac_c  = (__mT*) vv_fac_c.at(iswp).at(i_level);
            __mT* v_fac_f  = (__mT*) vv_fac_f.at(iswp).at(i_level);
            __mT* v_T0v    = (__mT*) vv_T0v.at(iswp).at(i_level);
            __mT* v_T0r    = (__mT*) vv_T0r.at(iswp).at(i_level);
            __mT* v_T0t    = (__mT*) vv_T0t.at(iswp).at(i_level);
            __mT* v_T0p    = (__mT*) vv_T0p.at(iswp).at(i_level);
            __mT* v_fun    = (__mT*) vv_fun.at(iswp).at(i_level);
            __mT* v_change = (__mT*) vv_change.at(iswp).at(i_level);

            // alias for dumped index
            int* dump_ijk   = vv_i__j__k__.at(iswp).at(i_level);
            int* dump_ip1jk = vv_ip1j__k__.at(iswp).at(i_level);
            int* dump_im1jk = vv_im1j__k__.at(iswp).at(i_level);
            int* dump_ijp1k = vv_i__jp1k__.at(iswp).at(i_level);
            int* dump_ijm1k = vv_i__jm1k__.at(iswp).at(i_level);
            int* dump_ijkp1 = vv_i__j__kp1.at(iswp).at(i_level);
            int* dump_ijkm1 = vv_i__j__km1.at(iswp).at(i_level);

            // load data of all nodes in one level on temporal aligned array
            for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {

                int i_vec = _i_vec * NSIMD;
                __mT v_c__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijk[i_vec]);
                __mT v_p__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ip1jk[i_vec]);
                __mT v_m__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_im1jk[i_vec]);
                __mT v__p_    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijp1k[i_vec]);
                __mT v__m_    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijm1k[i_vec]);
                __mT v___p    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkp1[i_vec]);
                __mT v___m    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkm1[i_vec]);

                // loop over all nodes in one level
                vect_stencil_1st_pre_simd(v_c__, \
                                          v_p__,    v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                          v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                          v_DP_inv, v_DT_inv, v_DR_inv, \
                                          v_DP_inv_half, v_DT_inv_half, v_DR_inv_half, \
                                          loc_I, loc_J, loc_K);

                // calculate updated value on c
                vect_stencil_1st_3rd_apre_simd(v_c__, v_fac_a[_i_vec], v_fac_b[_i_vec], v_fac_c[_i_vec], v_fac_f[_i_vec], \
                                               v_T0v[_i_vec], v_T0p[_i_vec]  , v_T0t[_i_vec]  , v_T0r[_i_vec]  , v_fun[_i_vec]  , v_change[_i_vec], \
                                               v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                               v_DP_inv, v_DT_inv, v_DR_inv);

                // store v_c__ to dump_c__
                _mmT_store_pT(dump_c__, v_c__);


                for (int i = 0; i < NSIMD; i++) {
                    if(i_vec+i>=n_nodes) break;

                    grid.tau_loc[dump_ijk[i_vec+i]] = dump_c__[i];
                }



            } // end of i_vec loop

            // mpi synchronization
            synchronize_all_sub();

        } // end of i_level loop

#elif USE_ARM_SVE

        svbool_t pg;
        //
        __mT v_DP_inv      = svdup_f64(1.0/dp);
        __mT v_DT_inv      = svdup_f64(1.0/dt);
        __mT v_DR_inv      = svdup_f64(1.0/dr);
        __mT v_DP_inv_half = svdup_f64(1.0/dp*0.5);
        __mT v_DT_inv_half = svdup_f64(1.0/dt*0.5);
        __mT v_DR_inv_half = svdup_f64(1.0/dr*0.5);

        // store stencil coefs
        __mT v_pp1;
        __mT v_pp2;
        __mT v_pt1;
        __mT v_pt2;
        __mT v_pr1;
        __mT v_pr2;

        __mT v_c__   ;
        __mT v_p__   ;
        __mT v_m__   ;
        __mT v__p_   ;
        __mT v__m_   ;
        __mT v___p   ;
        __mT v___m   ;

        __mT v_fac_a_ ;
        __mT v_fac_b_ ;
        __mT v_fac_c_ ;
        __mT v_fac_f_ ;
        __mT v_T0v_   ;
        __mT v_T0r_   ;
        __mT v_T0t_   ;
        __mT v_T0p_   ;
        __mT v_fun_   ;
        __mT v_change_;


        // measure time for only loop
        //auto start = std::chrono::high_resolution_clock::now();

        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc.at(i_level).size();
            //std::cout << "n_nodes = " << n_nodes << std::endl;

            int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

            // make alias to preloaded data
            CUSTOMREAL* v_fac_a  = vv_fac_a.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_b  = vv_fac_b.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_c  = vv_fac_c.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_f  = vv_fac_f.at(iswp).at(i_level);
            CUSTOMREAL* v_T0v    = vv_T0v.at(iswp).at(i_level);
            CUSTOMREAL* v_T0r    = vv_T0r.at(iswp).at(i_level);
            CUSTOMREAL* v_T0t    = vv_T0t.at(iswp).at(i_level);
            CUSTOMREAL* v_T0p    = vv_T0p.at(iswp).at(i_level);
            CUSTOMREAL* v_fun    = vv_fun.at(iswp).at(i_level);
            CUSTOMREAL* v_change = vv_change.at(iswp).at(i_level);

            // alias for dumped index
            uint64_t* dump_ijk   = vv_i__j__k__.at(iswp).at(i_level);
            uint64_t* dump_ip1jk = vv_ip1j__k__.at(iswp).at(i_level);
            uint64_t* dump_im1jk = vv_im1j__k__.at(iswp).at(i_level);
            uint64_t* dump_ijp1k = vv_i__jp1k__.at(iswp).at(i_level);
            uint64_t* dump_ijm1k = vv_i__jm1k__.at(iswp).at(i_level);
            uint64_t* dump_ijkp1 = vv_i__j__kp1.at(iswp).at(i_level);
            uint64_t* dump_ijkm1 = vv_i__j__km1.at(iswp).at(i_level);

            // load data of all nodes in one level on temporal aligned array
            for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {
                int i_vec = _i_vec * NSIMD;

                pg = svwhilelt_b64(i_vec, n_nodes);

                v_c__    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijk[i_vec]);
                v_p__    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ip1jk[i_vec]);
                v_m__    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_im1jk[i_vec]);
                v__p_    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijp1k[i_vec]);
                v__m_    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijm1k[i_vec]);
                v___p    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijkp1[i_vec]);
                v___m    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijkm1[i_vec]);

                // load v_iip, v_jjt, v_kkr
                v_fac_a_ = svld1_vnum_f64(pg, v_fac_a , _i_vec);
                v_fac_b_ = svld1_vnum_f64(pg, v_fac_b , _i_vec);
                v_fac_c_ = svld1_vnum_f64(pg, v_fac_c , _i_vec);
                v_fac_f_ = svld1_vnum_f64(pg, v_fac_f , _i_vec);
                v_T0v_   = svld1_vnum_f64(pg, v_T0v   , _i_vec);
                v_T0r_   = svld1_vnum_f64(pg, v_T0r   , _i_vec);
                v_T0t_   = svld1_vnum_f64(pg, v_T0t   , _i_vec);
                v_T0p_   = svld1_vnum_f64(pg, v_T0p   , _i_vec);
                v_fun_   = svld1_vnum_f64(pg, v_fun   , _i_vec);
                v_change_= svld1_vnum_f64(pg, v_change, _i_vec);

                // loop over all nodes in one level
                vect_stencil_1st_pre_simd(pg, \
                                          v_c__, \
                                          v_p__,    v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                          v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                          v_DP_inv, v_DT_inv, v_DR_inv, \
                                          v_DP_inv_half, v_DT_inv_half, v_DR_inv_half, \
                                          loc_I, loc_J, loc_K);

                // calculate updated value on c
                vect_stencil_1st_3rd_apre_simd(pg, v_c__, v_fac_a_, v_fac_b_, v_fac_c_, v_fac_f_, \
                                               v_T0v_, v_T0p_, v_T0t_, v_T0r_, v_fun_, v_change_, \
                                               v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                               v_DP_inv, v_DT_inv, v_DR_inv);

                // store v_c__ to dump_c__
                svst1_scatter_u64index_f64(pg, grid.tau_loc, svld1_u64(pg,&dump_ijk[i_vec]), v_c__);


            } // end of i_vec loop

            // mpi synchronization
            synchronize_all_sub();

        } // end of i_level loop

#endif // ifndef USE_SIMD

    } // end of if !use_gpu
    else { // if use_gpu

#if defined USE_CUDA

        // copy tau to device
        cuda_copy_tau_to_device(gpu_grid, grid.tau_loc);

        // run iteration
        cuda_run_iteration_forward(gpu_grid, iswp);

        // copy tau to host
        cuda_copy_tau_to_host(gpu_grid, grid.tau_loc);

#else // !defiend USE_CUDA
        // exit code
        std::cout << "Error: USE_CUDA is not defined" << std::endl;
        exit(1);
#endif

    } // end of if use_gpu

    // update boundary
    if (subdom_main) {
        calculate_boundary_nodes(grid);
    }
}


// ERROR index!!!!!!!!!!!!!
Iterator_level_3rd_order::Iterator_level_3rd_order(InputParams& IP, Grid& grid, Source& src, IO_utils& io, const std::string& src_name, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level(IP, grid, src, io, src_name, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_3rd_order::do_sweep(int iswp, Grid& grid, InputParams& IP){

    if(!use_gpu) {

#if !defined USE_SIMD

        // set sweep direction
        set_sweep_direction(iswp);

        int iip, jjt, kkr;
        int n_levels = ijk_for_this_subproc.size();

        for (int i_level = 0; i_level < n_levels; i_level++) {
            size_t n_nodes = ijk_for_this_subproc[i_level].size();

            for (size_t i_node = 0; i_node < n_nodes; i_node++) {

                V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

                if (r_dirc < 0) kkr = nr-kkr; //kk-1;
                else            kkr = kkr-1;  //nr-kk;
                if (t_dirc < 0) jjt = nt-jjt; //jj-1;
                else            jjt = jjt-1;  //nt-jj;
                if (p_dirc < 0) iip = np-iip; //ii-1;
                else            iip = iip-1;  //np-ii;

                //
                // calculate stencils
                //
                if (grid.is_changed[I2V(iip, jjt, kkr)]) {
                    calculate_stencil_3rd_order(grid, iip, jjt, kkr);
                } // is_changed == true
            } // end ijk

            // mpi synchronization
            synchronize_all_sub();

        } // end loop i_level

#elif USE_AVX512 || USE_AVX

        //
        __mT v_DP_inv      = _mmT_set1_pT(1.0/dp);
        __mT v_DT_inv      = _mmT_set1_pT(1.0/dt);
        __mT v_DR_inv      = _mmT_set1_pT(1.0/dr);
        __mT v_DP_inv_half = _mmT_set1_pT(1.0/dp*0.5);
        __mT v_DT_inv_half = _mmT_set1_pT(1.0/dt*0.5);
        __mT v_DR_inv_half = _mmT_set1_pT(1.0/dr*0.5);

        // store stencil coefs
        __mT v_pp1 = _mmT_set1_pT(0.0);
        __mT v_pp2 = _mmT_set1_pT(0.0);
        __mT v_pt1 = _mmT_set1_pT(0.0);
        __mT v_pt2 = _mmT_set1_pT(0.0);
        __mT v_pr1 = _mmT_set1_pT(0.0);
        __mT v_pr2 = _mmT_set1_pT(0.0);

        // measure time for only loop
        //auto start = std::chrono::high_resolution_clock::now();

        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc.at(i_level).size();
            //std::cout << "n_nodes = " << n_nodes << std::endl;

            int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

            // make alias to preloaded data
            __mT* v_iip    = (__mT*) vv_iip.at(iswp).at(i_level);
            __mT* v_jjt    = (__mT*) vv_jjt.at(iswp).at(i_level);
            __mT* v_kkr    = (__mT*) vv_kkr.at(iswp).at(i_level);

            __mT* v_fac_a  = (__mT*) vv_fac_a.at(iswp).at(i_level);
            __mT* v_fac_b  = (__mT*) vv_fac_b.at(iswp).at(i_level);
            __mT* v_fac_c  = (__mT*) vv_fac_c.at(iswp).at(i_level);
            __mT* v_fac_f  = (__mT*) vv_fac_f.at(iswp).at(i_level);
            __mT* v_T0v    = (__mT*) vv_T0v.at(iswp).at(i_level);
            __mT* v_T0r    = (__mT*) vv_T0r.at(iswp).at(i_level);
            __mT* v_T0t    = (__mT*) vv_T0t.at(iswp).at(i_level);
            __mT* v_T0p    = (__mT*) vv_T0p.at(iswp).at(i_level);
            __mT* v_fun    = (__mT*) vv_fun.at(iswp).at(i_level);
            __mT* v_change = (__mT*) vv_change.at(iswp).at(i_level);

            // alias for dumped index
            int* dump_ijk   = vv_i__j__k__.at(iswp).at(i_level);
            int* dump_ip1jk = vv_ip1j__k__.at(iswp).at(i_level);
            int* dump_im1jk = vv_im1j__k__.at(iswp).at(i_level);
            int* dump_ijp1k = vv_i__jp1k__.at(iswp).at(i_level);
            int* dump_ijm1k = vv_i__jm1k__.at(iswp).at(i_level);
            int* dump_ijkp1 = vv_i__j__kp1.at(iswp).at(i_level);
            int* dump_ijkm1 = vv_i__j__km1.at(iswp).at(i_level);
            int* dump_ip2jk = vv_ip2j__k__.at(iswp).at(i_level);
            int* dump_im2jk = vv_im2j__k__.at(iswp).at(i_level);
            int* dump_ijp2k = vv_i__jp2k__.at(iswp).at(i_level);
            int* dump_ijm2k = vv_i__jm2k__.at(iswp).at(i_level);
            int* dump_ijkp2 = vv_i__j__kp2.at(iswp).at(i_level);
            int* dump_ijkm2 = vv_i__j__km2.at(iswp).at(i_level); /////////////////////////////////////

            // load data of all nodes in one level on temporal aligned array
            for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {

                int i_vec = _i_vec * NSIMD;
                __mT v_c__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijk  [i_vec]);
                __mT v_p__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ip1jk[i_vec]);
                __mT v_m__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_im1jk[i_vec]);
                __mT v__p_    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijp1k[i_vec]);
                __mT v__m_    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijm1k[i_vec]);
                __mT v___p    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkp1[i_vec]);
                __mT v___m    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkm1[i_vec]);
                __mT v_pp____ = load_mem_gen_to_mTd(grid.tau_loc, &dump_ip2jk[i_vec]);
                __mT v_mm____ = load_mem_gen_to_mTd(grid.tau_loc, &dump_im2jk[i_vec]);
                __mT v___pp__ = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijp2k[i_vec]);
                __mT v___mm__ = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijm2k[i_vec]);
                __mT v_____pp = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkp2[i_vec]);
                __mT v_____mm = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkm2[i_vec]);

                // loop over all nodes in one level
                vect_stencil_3rd_pre_simd(v_iip[_i_vec], v_jjt[_i_vec], v_kkr[_i_vec], \
                                          v_c__, \
                                          v_p__,    v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                          v_pp____, v_mm____, v___pp__, v___mm__, v_____pp, v_____mm, \
                                          v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                          v_DP_inv, v_DT_inv, v_DR_inv, \
                                          v_DP_inv_half, v_DT_inv_half, v_DR_inv_half, \
                                          loc_I, loc_J, loc_K);

                // calculate updated value on c
                vect_stencil_1st_3rd_apre_simd(v_c__, v_fac_a[_i_vec], v_fac_b[_i_vec], v_fac_c[_i_vec], v_fac_f[_i_vec], \
                                               v_T0v[_i_vec], v_T0p[_i_vec]  , v_T0t[_i_vec]  , v_T0r[_i_vec]  , v_fun[_i_vec]  , v_change[_i_vec], \
                                               v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                               v_DP_inv, v_DT_inv, v_DR_inv);

                // store v_c__ to dump_c__
                _mmT_store_pT(dump_c__, v_c__);


                for (int i = 0; i < NSIMD; i++) {
                    if(i_vec+i>=n_nodes) break;

                    grid.tau_loc[dump_ijk[i_vec+i]] = dump_c__[i];
                }



            } // end of i_vec loop

            // mpi synchronization
            synchronize_all_sub(); // dead lock

        } // end of i_level loop

#elif USE_ARM_SVE
        //
        svbool_t pg;

        __mT v_DP_inv      = svdup_f64(1.0/dp);
        __mT v_DT_inv      = svdup_f64(1.0/dt);
        __mT v_DR_inv      = svdup_f64(1.0/dr);
        __mT v_DP_inv_half = svdup_f64(1.0/dp*0.5);
        __mT v_DT_inv_half = svdup_f64(1.0/dt*0.5);
        __mT v_DR_inv_half = svdup_f64(1.0/dr*0.5);

        // store stencil coefs
        __mT v_pp1;
        __mT v_pp2;
        __mT v_pt1;
        __mT v_pt2;
        __mT v_pr1;
        __mT v_pr2;

        __mT v_c__   ;
        __mT v_p__   ;
        __mT v_m__   ;
        __mT v__p_   ;
        __mT v__m_   ;
        __mT v___p   ;
        __mT v___m   ;
        __mT v_pp____;
        __mT v_mm____;
        __mT v___pp__;
        __mT v___mm__;
        __mT v_____pp;
        __mT v_____mm;

        __mT v_iip_   ;
        __mT v_jjt_   ;
        __mT v_kkr_   ;
        __mT v_fac_a_ ;
        __mT v_fac_b_ ;
        __mT v_fac_c_ ;
        __mT v_fac_f_ ;
        __mT v_T0v_   ;
        __mT v_T0r_   ;
        __mT v_T0t_   ;
        __mT v_T0p_   ;
        __mT v_fun_   ;
        __mT v_change_;


        // measure time for only loop
        //auto start = std::chrono::high_resolution_clock::now();

        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc.at(i_level).size();
            //std::cout << "n_nodes = " << n_nodes << std::endl;

            int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

            // make alias to preloaded data
            CUSTOMREAL* v_iip    = vv_iip.at(iswp).at(i_level);
            CUSTOMREAL* v_jjt    = vv_jjt.at(iswp).at(i_level);
            CUSTOMREAL* v_kkr    = vv_kkr.at(iswp).at(i_level);

            CUSTOMREAL* v_fac_a  = vv_fac_a.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_b  = vv_fac_b.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_c  = vv_fac_c.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_f  = vv_fac_f.at(iswp).at(i_level);
            CUSTOMREAL* v_T0v    = vv_T0v.at(iswp).at(i_level);
            CUSTOMREAL* v_T0r    = vv_T0r.at(iswp).at(i_level);
            CUSTOMREAL* v_T0t    = vv_T0t.at(iswp).at(i_level);
            CUSTOMREAL* v_T0p    = vv_T0p.at(iswp).at(i_level);
            CUSTOMREAL* v_fun    = vv_fun.at(iswp).at(i_level);
            CUSTOMREAL* v_change = vv_change.at(iswp).at(i_level);

            // alias for dumped index
            uint64_t* dump_ijk   = vv_i__j__k__.at(iswp).at(i_level);
            uint64_t* dump_ip1jk = vv_ip1j__k__.at(iswp).at(i_level);
            uint64_t* dump_im1jk = vv_im1j__k__.at(iswp).at(i_level);
            uint64_t* dump_ijp1k = vv_i__jp1k__.at(iswp).at(i_level);
            uint64_t* dump_ijm1k = vv_i__jm1k__.at(iswp).at(i_level);
            uint64_t* dump_ijkp1 = vv_i__j__kp1.at(iswp).at(i_level);
            uint64_t* dump_ijkm1 = vv_i__j__km1.at(iswp).at(i_level);
            uint64_t* dump_ip2jk = vv_ip2j__k__.at(iswp).at(i_level);
            uint64_t* dump_im2jk = vv_im2j__k__.at(iswp).at(i_level);
            uint64_t* dump_ijp2k = vv_i__jp2k__.at(iswp).at(i_level);
            uint64_t* dump_ijm2k = vv_i__jm2k__.at(iswp).at(i_level);
            uint64_t* dump_ijkp2 = vv_i__j__kp2.at(iswp).at(i_level);
            uint64_t* dump_ijkm2 = vv_i__j__km2.at(iswp).at(i_level);

            // load data of all nodes in one level on temporal aligned array
            for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {
                int i_vec = _i_vec * NSIMD;

                pg = svwhilelt_b64(i_vec, n_nodes);

                v_c__    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijk  [i_vec]);
                v_p__    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ip1jk[i_vec]);
                v_m__    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_im1jk[i_vec]);
                v__p_    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijp1k[i_vec]);
                v__m_    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijm1k[i_vec]);
                v___p    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijkp1[i_vec]);
                v___m    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijkm1[i_vec]);
                v_pp____ = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ip2jk[i_vec]);
                v_mm____ = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_im2jk[i_vec]);
                v___pp__ = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijp2k[i_vec]);
                v___mm__ = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijm2k[i_vec]);
                v_____pp = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijkp2[i_vec]);
                v_____mm = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijkm2[i_vec]);

                // load v_iip, v_jjt, v_kkr
                v_iip_   = svld1_vnum_f64(pg, v_iip   , _i_vec);
                v_jjt_   = svld1_vnum_f64(pg, v_jjt   , _i_vec);
                v_kkr_   = svld1_vnum_f64(pg, v_kkr   , _i_vec);
                v_fac_a_ = svld1_vnum_f64(pg, v_fac_a , _i_vec);
                v_fac_b_ = svld1_vnum_f64(pg, v_fac_b , _i_vec);
                v_fac_c_ = svld1_vnum_f64(pg, v_fac_c , _i_vec);
                v_fac_f_ = svld1_vnum_f64(pg, v_fac_f , _i_vec);
                v_T0v_   = svld1_vnum_f64(pg, v_T0v   , _i_vec);
                v_T0r_   = svld1_vnum_f64(pg, v_T0r   , _i_vec);
                v_T0t_   = svld1_vnum_f64(pg, v_T0t   , _i_vec);
                v_T0p_   = svld1_vnum_f64(pg, v_T0p   , _i_vec);
                v_fun_   = svld1_vnum_f64(pg, v_fun   , _i_vec);
                v_change_= svld1_vnum_f64(pg, v_change, _i_vec);

                // loop over all nodes in one level
                vect_stencil_3rd_pre_simd(pg, v_iip_, v_jjt_, v_kkr_, \
                                          v_c__, \
                                          v_p__,    v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                          v_pp____, v_mm____, v___pp__, v___mm__, v_____pp, v_____mm, \
                                          v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                          v_DP_inv, v_DT_inv, v_DR_inv, \
                                          v_DP_inv_half, v_DT_inv_half, v_DR_inv_half, \
                                          loc_I, loc_J, loc_K);

                //// calculate updated value on c
                vect_stencil_1st_3rd_apre_simd(pg, v_c__, v_fac_a_, v_fac_b_, v_fac_c_, v_fac_f_, \
                                               v_T0v_, v_T0p_, v_T0t_, v_T0r_, v_fun_, v_change_, \
                                               v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                               v_DP_inv, v_DT_inv, v_DR_inv);

                // store v_c__ to dump_c__
                svst1_scatter_u64index_f64(pg, grid.tau_loc, svld1_u64(pg,&dump_ijk[i_vec]), v_c__);
            } // end of i_vec loop

            // mpi synchronization
            synchronize_all_sub();

        } // end of i_level loop

#endif // ifndef USE_SIMD

    } // end of if !use_gpu
    else { // if use_gpu

#if defined USE_CUDA

        // copy tau to device
        cuda_copy_tau_to_device(gpu_grid, grid.tau_loc);

        // run iteration
        cuda_run_iteration_forward(gpu_grid, iswp);

        // copy tau to host
        cuda_copy_tau_to_host(gpu_grid, grid.tau_loc);

#else // !defiend USE_CUDA
        // exit code
        std::cout << "Error: USE_CUDA is not defined" << std::endl;
        exit(1);
#endif

    } // end of if use_gpu

    // update boundary
    if (subdom_main) {
        calculate_boundary_nodes(grid);
    }

}


Iterator_level_1st_order_upwind::Iterator_level_1st_order_upwind(InputParams& IP, Grid& grid, Source& src, IO_utils& io, const std::string& src_name, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level(IP, grid, src, io, src_name, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}

void Iterator_level_1st_order_upwind::do_sweep(int iswp, Grid& grid, InputParams& IP){

    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-1-kkr; 
            else            kkr = kkr;
            if (t_dirc < 0) jjt = nt-1-jjt; 
            else            jjt = jjt;  
            if (p_dirc < 0) iip = np-1-iip; 
            else            iip = iip;  

            //
            // calculate stencils
            //
            if (grid.is_changed[I2V(iip, jjt, kkr)]) {
                calculate_stencil_1st_order_upwind(grid, iip, jjt, kkr);
            } // is_changed == true
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level
}

// ERROR index!!!!!!!!!!!!!
Iterator_level_1st_order_tele::Iterator_level_1st_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, const std::string& src_name, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level_tele(IP, grid, src, io, src_name, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}

void Iterator_level_1st_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){

    if(!use_gpu) {

//#if !defined USE_SIMD

        // set sweep direction
        set_sweep_direction(iswp);

        int iip, jjt, kkr;
        int n_levels = ijk_for_this_subproc.size();

        for (int i_level = 0; i_level < n_levels; i_level++) {
            size_t n_nodes = ijk_for_this_subproc[i_level].size();

            for (size_t i_node = 0; i_node < n_nodes; i_node++) {

                V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

                if (r_dirc < 0) kkr = nr-kkr-1;
                //else            kkr = kkr;
                if (t_dirc < 0) jjt = nt-jjt-1;
                //else            jjt = jjt;
                if (p_dirc < 0) iip = np-iip-1;
                //else            iip = iip;

                //
                // calculate stencils
                //
                if (iip != 0 && iip != np-1 && jjt != 0 && jjt != nt-1 && kkr != 0 && kkr != nr-1) {
                    // calculate stencils
                    calculate_stencil_1st_order_tele(grid, iip, jjt, kkr);
                } else {
                    // update boundary
                    calculate_boundary_nodes_tele(grid, iip, jjt, kkr);
                }
            } // end ijk

            // mpi synchronization
            synchronize_all_sub();

        } // end loop i_level

/*
#elif USE_AVX512 || USE_AVX

        // preload constants
        __mT v_DP_inv      = _mmT_set1_pT(1.0/dp);
        __mT v_DT_inv      = _mmT_set1_pT(1.0/dt);
        __mT v_DR_inv      = _mmT_set1_pT(1.0/dr);
        __mT v_DP_inv_half = _mmT_set1_pT(1.0/dp*0.5);
        __mT v_DT_inv_half = _mmT_set1_pT(1.0/dt*0.5);
        __mT v_DR_inv_half = _mmT_set1_pT(1.0/dr*0.5);

        // store stencil coefs
        __mT v_pp1;
        __mT v_pp2;
        __mT v_pt1;
        __mT v_pt2;
        __mT v_pr1;
        __mT v_pr2;

        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc.at(i_level).size();

            int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

            // make alias to preloaded data
            __mT* v_iip    = (__mT*) vv_iip.at(iswp).at(i_level);
            __mT* v_jjt    = (__mT*) vv_jjt.at(iswp).at(i_level);
            __mT* v_kkr    = (__mT*) vv_kkr.at(iswp).at(i_level);

            __mT* v_fac_a  = (__mT*) vv_fac_a.at(iswp).at(i_level);
            __mT* v_fac_b  = (__mT*) vv_fac_b.at(iswp).at(i_level);
            __mT* v_fac_c  = (__mT*) vv_fac_c.at(iswp).at(i_level);
            __mT* v_fac_f  = (__mT*) vv_fac_f.at(iswp).at(i_level);
            __mT* v_fun    = (__mT*) vv_fun.at(iswp).at(i_level);
            __mT* v_change = (__mT*) vv_change.at(iswp).at(i_level);

            // alias for dumped index
            int* dump_ijk   = vv_i__j__k__.at(iswp).at(i_level);
            int* dump_ip1jk = vv_ip1j__k__.at(iswp).at(i_level);
            int* dump_im1jk = vv_im1j__k__.at(iswp).at(i_level);
            int* dump_ijp1k = vv_i__jp1k__.at(iswp).at(i_level);
            int* dump_ijm1k = vv_i__jm1k__.at(iswp).at(i_level);
            int* dump_ijkp1 = vv_i__j__kp1.at(iswp).at(i_level);
            int* dump_ijkm1 = vv_i__j__km1.at(iswp).at(i_level);
            int* dump_ip2jk = vv_ip2j__k__.at(iswp).at(i_level);
            int* dump_im2jk = vv_im2j__k__.at(iswp).at(i_level);
            int* dump_ijp2k = vv_i__jp2k__.at(iswp).at(i_level);
            int* dump_ijm2k = vv_i__jm2k__.at(iswp).at(i_level);
            int* dump_ijkp2 = vv_i__j__kp2.at(iswp).at(i_level);
            int* dump_ijkm2 = vv_i__j__km2.at(iswp).at(i_level);

            // load data of all nodes in one level on temporal aligned array
            for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {

                int i_vec = _i_vec * NSIMD;
                __mT v_c__    = load_mem_gen_to_mTd(grid.T_loc,   &dump_ijk[i_vec]);
                __mT v_p__    = load_mem_gen_to_mTd(grid.T_loc,   &dump_ip1jk[i_vec]);
                __mT v_m__    = load_mem_gen_to_mTd(grid.T_loc,   &dump_im1jk[i_vec]);
                __mT v__p_    = load_mem_gen_to_mTd(grid.T_loc,   &dump_ijp1k[i_vec]);
                __mT v__m_    = load_mem_gen_to_mTd(grid.T_loc,   &dump_ijm1k[i_vec]);
                __mT v___p    = load_mem_gen_to_mTd(grid.T_loc,   &dump_ijkp1[i_vec]);
                __mT v___m    = load_mem_gen_to_mTd(grid.T_loc,   &dump_ijkm1[i_vec]);
                __mT v_pp____ = load_mem_gen_to_mTd(grid.T_loc, &dump_ip2jk[i_vec]);
                __mT v_mm____ = load_mem_gen_to_mTd(grid.T_loc, &dump_im2jk[i_vec]);
                __mT v___pp__ = load_mem_gen_to_mTd(grid.T_loc, &dump_ijp2k[i_vec]);
                __mT v___mm__ = load_mem_gen_to_mTd(grid.T_loc, &dump_ijm2k[i_vec]);
                __mT v_____pp = load_mem_gen_to_mTd(grid.T_loc, &dump_ijkp2[i_vec]);
                __mT v_____mm = load_mem_gen_to_mTd(grid.T_loc, &dump_ijkm2[i_vec]);


                // loop over all nodes in one level (this routine for teleseismic is same with that for local)
                vect_stencil_1st_pre_simd(v_iip[_i_vec], v_jjt[_i_vec], v_kkr[_i_vec], \
                                          v_c__, \
                                          v_p__,    v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                          v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                          v_DP_inv, v_DT_inv, v_DR_inv, \
                                          v_DP_inv_half, v_DT_inv_half, v_DR_inv_half, \
                                          loc_I, loc_J, loc_K);

                // calculate updated value on c
                vect_stencil_1st_3rd_apre_simd_tele(v_c__, v_fac_a[_i_vec], v_fac_b[_i_vec], v_fac_c[_i_vec], v_fac_f[_i_vec], \
                                               v_fun[_i_vec], v_change[_i_vec], \
                                               v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                               v_DP_inv, v_DT_inv, v_DR_inv);

                // calculate the values on boudaries
                calculate_boundary_nodes_tele_simd(v_iip[_i_vec], v_jjt[_i_vec], v_kkr[_i_vec], \
                                                   v_c__, \
                                                   v_p__,   v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                                   v_pp____, v_mm____, v___pp__, v___mm__, v_____pp, v_____mm, \
                                                   v_change[_i_vec], \
                                                   loc_I, loc_J, loc_K);

                // store v_c__ to dump_c__
                _mmT_store_pT(dump_c__, v_c__);

                for (int i = 0; i < NSIMD; i++) {
                    if(i_vec+i>=n_nodes) break;

                    grid.tau_loc[dump_ijk[i_vec+i]] = dump_c__[i];
                }



            } // end of i_vec loop

            // mpi synchronization
            synchronize_all_sub();

        } // end of i_level loop



#elif USE_ARM_SVE

#endif // ifndef USE_SIMD
*/

    } // end of if !use_gpu
    else { // if use_gpu

#if defined USE_CUDA

        //// copy tau to device
        //cuda_copy_tau_to_device(gpu_grid, grid.tau_loc);

        //// run iteration
        //cuda_run_iteration_forward_tele(gpu_grid, iswp);

        //// copy tau to host
        //cuda_copy_tau_to_host(gpu_grid, grid.tau_loc);

#else // !defiend USE_CUDA
        // exit code
        std::cout << "Error: USE_CUDA is not defined" << std::endl;
        exit(1);
#endif

    } // end of if use_gpu


}

// ERROR index!!!!!!!!!!!!!
Iterator_level_3rd_order_tele::Iterator_level_3rd_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, const std::string& src_name, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level_tele(IP, grid, src, io, src_name, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_3rd_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){

    if(!use_gpu) {

//#if !defined USE_SIMD

        // set sweep direction
        set_sweep_direction(iswp);

        int iip, jjt, kkr;
        int n_levels = ijk_for_this_subproc.size();

        for (int i_level = 0; i_level < n_levels; i_level++) {
            size_t n_nodes = ijk_for_this_subproc[i_level].size();

            for (size_t i_node = 0; i_node < n_nodes; i_node++) {

                V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

                if (r_dirc < 0) kkr = nr-kkr-1;
                //else            kkr = kkr;
                if (t_dirc < 0) jjt = nt-jjt-1;
                //else            jjt = jjt;
                if (p_dirc < 0) iip = np-iip-1;
                //else            iip = iip;

                //
                // calculate stencils
                //
                if (iip != 0 && iip != np-1 && jjt != 0 && jjt != nt-1 && kkr != 0 && kkr != nr-1) {
                    // calculate stencils
                    calculate_stencil_3rd_order_tele(grid, iip, jjt, kkr);
                } else {
                    // update boundary
                    calculate_boundary_nodes_tele(grid, iip, jjt, kkr);
                }
            } // end ijk

            // mpi synchronization
            synchronize_all_sub(); ///////////////////////////////////////

        } // end loop i_level

//#elif USE_AVX512 || USE_AVX
//
//#elif USE_ARM_SVE
//
//#endif // ifndef USE_SIMD

    } // end of if !use_gpu
    else { // if use_gpu

#if defined USE_CUDA

        //// copy tau to device
        //cuda_copy_tau_to_device(gpu_grid, grid.tau_loc);

        //// run iteration
        //cuda_run_iteration_forward_tele(gpu_grid, iswp);

        //// copy tau to host
        //cuda_copy_tau_to_host(gpu_grid, grid.tau_loc);

#else // !defiend USE_CUDA
        // exit code
        std::cout << "Error: USE_CUDA is not defined" << std::endl;
        exit(1);
#endif

    } // end of if use_gpu


}


Iterator_level_1st_order_upwind_tele::Iterator_level_1st_order_upwind_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, const std::string& src_name, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level_tele(IP, grid, src, io, src_name, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}

void Iterator_level_1st_order_upwind_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){

    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-1-kkr; 
            else            kkr = kkr;
            if (t_dirc < 0) jjt = nt-1-jjt; 
            else            jjt = jjt;  
            if (p_dirc < 0) iip = np-1-iip; 
            else            iip = iip;  

            //
            // calculate stencils
            //
            // if (iip != 0 && iip != np-1 && jjt != 0 && jjt != nt-1 && kkr != 0) {   // top layer is not fixed, otherwise, the top layer will be 2000
            //     calculate_stencil_1st_order_upwind_tele(grid, iip, jjt, kkr);       // no need to consider the boundary for upwind scheme
            // }

            calculate_stencil_1st_order_upwind_tele(grid, iip, jjt, kkr);
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level
}
