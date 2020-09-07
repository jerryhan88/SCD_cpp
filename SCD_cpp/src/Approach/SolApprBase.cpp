//
//  SolApprBase.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 13/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/SolApprBase.hpp"


void BaseMM::generate_y_inv(char vType) {
    if (lp_algo == "DSM") {
        y_ary = new_ak_inv(prob, env, vType);
    } else if (lp_algo == "IPM"){
        new_ak_inv(prob, env, vType, y_map);
    } else {
        assert(false);
    }
}

void BaseMM::generate_z_inv(char vType) {
    if (lp_algo == "DSM") {
        z_ary = new_aek_inv(prob, env, vType);
    } else if (lp_algo == "IPM"){
        new_aek_inv(prob, env, vType, z_map);
    } else {
        assert(false);
    }
}

void BaseMM::generate_x_inv(char vType) {
    if (lp_algo == "DSM") {
        x_ary = new_aeij_inv(prob, env, vType, bool_x_aeij);
    } else if (lp_algo == "IPM"){
        new_aeij_inv(prob, env, vType, bool_x_aeij, x_map);
    } else {
        assert(false);
    }
}

void BaseMM::generate_u_inv() {
    if (lp_algo == "DSM") {
        u_ary = new_aei_inv(prob, env);
    } else if (lp_algo == "IPM"){
        new_aei_inv(prob, env, u_map);
    } else {
        assert(false);
    }
}


void BaseMM::delete_y_inv() {
    if (lp_algo == "DSM") {
        del_ak_inv(prob, y_ary);
    } else if (lp_algo == "IPM"){
        y_map.clear();
    }
}

void BaseMM::delete_z_inv() {
    if (lp_algo == "DSM") {
        del_aek_inv(prob, z_ary);
    } else if (lp_algo == "IPM"){
        z_map.clear();
    }
}

void BaseMM::delete_x_inv() {
    if (lp_algo == "DSM") {
        del_aeij_inv(prob, x_ary);
    } else if (lp_algo == "IPM"){
        x_map.clear();
    }
}

void BaseMM::delete_u_inv() {
    if (lp_algo == "DSM") {
        del_aei_inv(prob, u_ary);
    } else if (lp_algo == "IPM"){
        u_map.clear();
    }
}


IloNumVar BaseMM::get_y_inv(int a, int k) {
    if (lp_algo == "DSM") {
        return y_ary[a][k];
    } else {
        return y_map[std::make_tuple(a, k)];
    }
}

IloNumVar BaseMM::get_z_inv(int a, int e, int k) {
    if (lp_algo == "DSM") {
        return z_ary[a][e][k];
    } else {
        return z_map[std::make_tuple(a, e, k)];
    }
}

IloNumVar BaseMM::get_x_inv(int a, int e, int i, int j) {
    if (lp_algo == "DSM") {
        return x_ary[a][e][i][j];
    } else {
        return x_map[std::make_tuple(a, e, i, j)];
    }
}

IloNumVar BaseMM::get_u_inv(int a, int e, int i) {
    if (lp_algo == "DSM") {
        return u_ary[a][e][i];
    } else {
        return u_map[std::make_tuple(a, e, i)];
    }
}


void BaseMM::build_baseModel() {
    def_objF();
    def_ETA_cnsts();
    def_RUT_cnsts();
    def_COM_cnsts();
}


void BaseMM::preprocessing() {
    int a, e, i, j, o, d;
    std::set<iiiiTup> invalid_dvX;
    //
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i0 = 0; i0 < prob->R_ae[a][e].size() - 1; i0++) {
                int n0 = prob->R_ae[a][e][i0];
                double al_n0 = prob->al_i[n0];
                int n2 = prob->R_ae[a][e][i0 + 1];
                double be_n2 = prob->be_i[n2];
                for (int n1: prob->P_ae[a][e]) {
                    if (al_n0 + prob->t_ij[n0][n1] + prob->t_ij[n1][n2] > be_n2) {
                        invalid_dvX.insert(std::make_tuple(a, e, n0, n1));
                        invalid_dvX.insert(std::make_tuple(a, e, n1, n2));
                    }
                }
                for (int n1: prob->D_ae[a][e]) {
                    if (al_n0 + prob->t_ij[n0][n1] + prob->t_ij[n1][n2] > be_n2) {
                        invalid_dvX.insert(std::make_tuple(a, e, n0, n1));
                        invalid_dvX.insert(std::make_tuple(a, e, n1, n2));
                    }
                }
            }
        }
    }
    //
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    if (i == prob->d_ae[a][e] && j == prob->o_ae[a][e]) {
                        continue;
                    }
                    if (prob->al_i[i] + prob->t_ij[i][j] > prob->be_i[j]) {
                        invalid_dvX.insert(std::make_tuple(a, e, i, j));
                    }
                }
            }
        }
    }
    //
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            o = prob->o_ae[a][e];
            for (int i: prob->N_ae[a][e]) {
                if (i == prob->d_ae[a][e]) {
                    continue;
                }
                invalid_dvX.insert(std::make_tuple(a, e, i, o));
            }
            d = prob->d_ae[a][e];
            for (int j: prob->N_ae[a][e]) {
                if (j == prob->o_ae[a][e]) {
                    continue;
                }
                invalid_dvX.insert(std::make_tuple(a, e, d, j));
            }
        }
    }
    //
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k : prob->K_ae[a][e]) {
                int i = prob->n_k[k];
                int j = prob->h_k[k];
                invalid_dvX.insert(std::make_tuple(a, e, i, j));
            }
        }
    }
    //
    for (iiiiTup ele: invalid_dvX) {
        std::tie(a, e, i, j) = ele;
        bool_x_aeij[a][e][i][i] = false;
    }
}

void BaseMM::def_objF() {
    IloExpr objF(env);
    for (int k: prob->K) {
        for (int a : prob->A) {
            objF += prob->r_k[k] * get_y_inv(a, k);
            for (int e: prob->E_a[a]) {
                objF -= prob->r_k[k] * prob->p_ae[a][e] * get_z_inv(a, e, k);
            }
        }
    }
    baseModel->add(IloMaximize(env, objF));
    objF.end();
}

void BaseMM::def_ETA_cnsts() {
    //
    // Evaluation of the Task Assignment
    //
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    for (int k : prob->K) {
        linExpr.clear();
        sprintf(buf, "TAS(%d)", k);  // Task Assignment
        for (int a : prob->A) {
            linExpr += get_y_inv(a, k);
        }
        cnsts.add(linExpr <= 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                if (prob->K_ae[a][e].find(k) == prob->K_ae[a][e].end()) {
                    linExpr.clear();
                    sprintf(buf, "ITA(%d)(%d)(%d)", a, e, k);  // Infeasible Tasks Assignment
                    linExpr += get_y_inv(a, k);
                    cnsts.add(linExpr == 0);
                    cnsts[cnsts.getSize() - 1].setName(buf);
                }
            }
        }
    }
    //
    for (int a: prob->A) {
        linExpr.clear();
        sprintf(buf, "VL(%d)", a);  // Volume Limit
        for (int k: prob->K) {
            linExpr += prob->v_k[k] * get_y_inv(a, k);
        }
        cnsts.add(linExpr <= prob->v_a[a]);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    for (int a: prob->A) {
        linExpr.clear();
        sprintf(buf, "WL(%d)", a);  // Weight Limit
        for (int k: prob->K) {
            linExpr += prob->w_k[k] * get_y_inv(a, k);
        }
        cnsts.add(linExpr <= prob->w_a[a]);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                linExpr.clear();
                sprintf(buf, "TAC(%d)(%d)(%d)", a, e, k);  // Task Accomplishment
                linExpr += get_z_inv(a, e, k);
                linExpr -= get_y_inv(a, k);
                cnsts.add(linExpr <= 0);
                cnsts[cnsts.getSize() - 1].setName(buf);
            }
        }
    }
    baseModel->add(cnsts);
}

void BaseMM::def_RUT_cnsts() {
    // Routing constraints
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            def_FC_cnsts_aeGiven(a, e);
            def_AT_cnsts_aeGiven(a, e);
        }
    }
}

void BaseMM::def_FC_cnsts_aeGiven(int a, int e) {
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    int i, j, o, d;
    //
    // Initiate flow
    //
    linExpr.clear();
    sprintf(buf, "iFO(%d)(%d)", a, e);
    o = prob->o_ae[a][e];
    for (int j0: prob->N_ae[a][e]) {
        if (!bool_x_aeij[a][e][o][j0]) {
            continue;
        }
        linExpr += get_x_inv(a, e, o, j0);
    }
    cnsts.add(linExpr == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "iFD(%d)(%d)", a, e);
    d = prob->d_ae[a][e];
    for (int i0: prob->N_ae[a][e]) {
        if (!bool_x_aeij[a][e][i0][d]) {
            continue;
        }
        linExpr += get_x_inv(a, e, i0, d);
    }
    cnsts.add(linExpr == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    for (int i: prob->R_ae[a][e]) {
        if (i == prob->o_ae[a][e] || i == prob->d_ae[a][e]) {
            continue;
        }
        linExpr.clear();
        sprintf(buf, "iFR1(%d)(%d)(%d)", a, e, i);
        for (int j: prob->N_ae[a][e]) {
            if (j == i) continue;
            linExpr += get_x_inv(a, e, i, j);
        }
        cnsts.add(linExpr == 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
        //
        linExpr.clear();
        sprintf(buf, "iFR2(%d)(%d)(%d)", a, e, i);
        for (int j: prob->N_ae[a][e]) {
            if (j == i) continue;
            linExpr += get_x_inv(a, e, j, i);
        }
        cnsts.add(linExpr == 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    sprintf(buf, "CF(%d)(%d)", a, e);  // Circular Flow for the branch-and-cut algorithm
    cnsts.add(get_x_inv(a, e, prob->d_ae[a][e], prob->o_ae[a][e]) == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    for (int i0: prob->N_ae[a][e]) {
        if (!bool_x_aeij[a][e][i0][i0]) {
            continue;
        }
        sprintf(buf, "NS(%d)(%d)(%d)", a, e, i0);  // No Self Flow; tightening bounds
        cnsts.add(get_x_inv(a, e, i0, i0) == 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    // Flow about delivery nodes; only when the warehouse visited
    //
    for (int k: prob->K_ae[a][e]) {
        linExpr.clear();
        sprintf(buf, "tFC(%d)(%d)(%d)", a, e, k);
        i = prob->n_k[k];
        for (int j0: prob->N_ae[a][e]) {
            if (!bool_x_aeij[a][e][i][j0]) {
                continue;
            }
            linExpr += get_x_inv(a, e, i, j0);
        }
        j = prob->h_k[k];
        for (int i0: prob->N_ae[a][e]) {
            if (!bool_x_aeij[a][e][i0][j]) {
                continue;
            }
            linExpr -= get_x_inv(a, e, i0, j);
        }
        cnsts.add(linExpr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    // Flow conservation
    //
    for (int i0: prob->PD_ae[a][e]) {
        linExpr.clear();
        sprintf(buf, "FC_1(%d)(%d)(%d)", a, e, i0);
        for (int j0: prob->N_ae[a][e]) {
            if (!bool_x_aeij[a][e][i0][j0]) {
                continue;
            }
            linExpr += get_x_inv(a, e, i0, j0);
        }
        cnsts.add(linExpr <= 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
        //
        sprintf(buf, "FC(%d)(%d)(%d)", a, e, i0);
        for (int j0: prob->N_ae[a][e]) {
            if (!bool_x_aeij[a][e][j0][i0]) {
                continue;
            }
            linExpr -= get_x_inv(a, e, j0, i0);
        }
        cnsts.add(linExpr == 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    linExpr.end();
    baseModel->add(cnsts);
}

void BaseMM::def_AT_cnsts_aeGiven(int a, int e) {
    //
    // Constraints associated with the arrival time
    //
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    // Initialize the arrival time on the origin node
    //
    sprintf(buf, "IA(%d)(%d)", a, e);
    cnsts.add(get_u_inv(a, e, prob->o_ae[a][e]) == prob->al_i[prob->o_ae[a][e]]);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    // Arrival time calculation
    //
    for (int i: prob->N_ae[a][e]) {
        for (int j: prob->N_ae[a][e]) {
            if (i == prob->d_ae[a][e] && j == prob->o_ae[a][e]) {
                continue;
            }
            if (!bool_x_aeij[a][e][i][j]) {
                continue;
            }
            linExpr.clear();
            sprintf(buf, "AT(%d)(%d)(%d)(%d)", a, e, i, j);
            linExpr += get_u_inv(a, e, i) + prob->t_ij[i][j];
            linExpr -= get_u_inv(a, e, j) + prob->M * (1 - get_x_inv(a, e, i, j));
            cnsts.add(linExpr <= 0);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
    }
    //
    // Time Window
    //
    for (int i: prob->N_ae[a][e]) {
        sprintf(buf, "TW1(%d)(%d)(%d)", a, e, i);
        cnsts.add(prob->al_i[i] <= get_u_inv(a, e, i));
        cnsts[cnsts.getSize() - 1].setName(buf);
        //
        sprintf(buf, "TW2(%d)(%d)(%d)", a, e, i);
        cnsts.add(get_u_inv(a, e, i) <= prob->be_i[i]);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    // Routine route preservation
    //
    for (int i: prob->R_ae[a][e]) {
        for (int j: prob->R_ae[a][e]) {
            linExpr.clear();
            sprintf(buf, "RR_P(%d)(%d)(%d)(%d)", a, e, i, j);
            linExpr += prob->c_aeij[a][e][i][j] * get_u_inv(a, e, i);
            linExpr -= get_u_inv(a, e, j);
            cnsts.add(linExpr <= 0);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
    }
    //
    // Warehouse and Delivery Sequence
    //
    for (int k: prob->K_ae[a][e]) {
        linExpr.clear();
        sprintf(buf, "WD_S(%d)(%d)(%d)", a, e, k);
        linExpr += get_u_inv(a, e, prob->h_k[k]);
        linExpr -= get_u_inv(a, e, prob->n_k[k]);
        linExpr -= prob->M;
        for (int j: prob->N_ae[a][e]) {
            if (!bool_x_aeij[a][e][prob->n_k[k]][j]) {
                continue;
            }
            linExpr += prob->M * get_x_inv(a, e, prob->n_k[k], j);
        }
        cnsts.add(linExpr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    // Volume Limit
    //
    linExpr.clear();
    sprintf(buf, "VL(%d)(%d)", a, e);
    for (int k: prob->K_ae[a][e]) {
        for (int j: prob->N_ae[a][e]) {
            if (!bool_x_aeij[a][e][j][prob->n_k[k]]) {
                continue;
            }
            linExpr += prob->v_k[k] * get_x_inv(a, e, j, prob->n_k[k]);
        }
    }
    cnsts.add(linExpr <= prob->v_a[a]);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    // Weight Limit
    //
    linExpr.clear();
    sprintf(buf, "WL(%d)(%d)", a, e);
    for (int k: prob->K_ae[a][e]) {
        for (int j: prob->N_ae[a][e]) {
            if (!bool_x_aeij[a][e][j][prob->n_k[k]]) {
                continue;
            }
            linExpr += prob->w_k[k] * get_x_inv(a, e, j, prob->n_k[k]);
        }
    }
    cnsts.add(linExpr <= prob->w_a[a]);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    // Time Limit
    //
    linExpr.clear();
    sprintf(buf, "TL(%d)(%d)", a, e);
    linExpr.clear();
    for (int i: prob->N_ae[a][e]) {
        for (int j: prob->N_ae[a][e]) {
            if (i == prob->d_ae[a][e] and j == prob->o_ae[a][e]) {
                continue;
            }
            if (!bool_x_aeij[a][e][i][j]) {
                continue;
            }
            linExpr += prob->t_ij[i][j] * get_x_inv(a, e, i, j);
        }
    }
    cnsts.add(linExpr <= prob->u_ae[a][e]);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.end();
    baseModel->add(cnsts);
}

void BaseMM::def_COM_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
    IloExpr linExpr(env);
    //
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                linExpr.clear();
                sprintf(buf, "CC(%d)(%d)(%d)", a, e, k);
                linExpr += get_y_inv(a, k);
                for (int j: prob->N_ae[a][e]) {
                    if (!bool_x_aeij[a][e][j][prob->n_k[k]]) {
                        continue;
                    }
                    linExpr -= get_x_inv(a, e, j, prob->n_k[k]);
                }
                linExpr -= get_z_inv(a, e, k);
                COM_cnsts->add(linExpr <= 0);
                (*COM_cnsts)[COM_cnsts->getSize() - 1].setName(buf);
                COM_cnsts_index[a][e][k] = COM_cnsts->getSize() - 1;
            }
        }
    }
    //
    linExpr.end();
    baseModel->add(*COM_cnsts);
}
