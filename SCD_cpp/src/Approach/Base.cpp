//
//  Base.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 17/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Base.hpp"

void Base::build_baseModel() {
    def_objF();
    def_ETA_cnsts(prob, env, y_ak, z_aek, baseModel);
    def_RUT_cnsts();
    def_COM_cnsts();
}

void Base::def_objF() {
    IloExpr objF(env);
    for (int k: prob->K) {
        for (int a : prob->A) {
            objF += prob->r_k[k] * y_ak[a][k];
            for (int e: prob->E_a[a]) {
                objF -= prob->r_k[k] * prob->p_ae[a][e] * z_aek[a][e][k];
            }
        }
    }
    baseModel->add(IloMaximize(env, objF));
    objF.end();
}

void Base::def_RUT_cnsts() {
    // Routing constraints
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            def_FC_cnsts_aeGiven(a, e);
            def_AT_cnsts_aeGiven(a, e);
        }
    }
}

void Base::def_FC_cnsts_aeGiven(int a, int e) {
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    // Initiate flow
    //
    linExpr.clear();
    sprintf(buf, "iFO(%d)(%d)", a, e);
    for (int j: prob->N_ae[a][e]) {
        linExpr += x_aeij[a][e][prob->o_ae[a][e]][j];
    }
    cnsts.add(linExpr == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "iFD(%d)(%d)", a, e);
    for (int j: prob->N_ae[a][e]) {
        linExpr += x_aeij[a][e][j][prob->d_ae[a][e]];
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
            linExpr += x_aeij[a][e][i][j];
        }
        cnsts.add(linExpr == 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
        //
        linExpr.clear();
        sprintf(buf, "iFR2(%d)(%d)(%d)", a, e, i);
        for (int j: prob->N_ae[a][e]) {
            if (j == i) continue;
            linExpr += x_aeij[a][e][j][i];
        }
        cnsts.add(linExpr == 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    sprintf(buf, "CF(%d)(%d)", a, e);  // Circular Flow for the branch-and-cut algorithm
    cnsts.add(x_aeij[a][e][prob->d_ae[a][e]][prob->o_ae[a][e]] == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "NS(%d)(%d)", a, e);  // No Self Flow; tightening bounds
    for (int i: prob->N_ae[a][e]) {
        linExpr += x_aeij[a][e][i][i];
    }
    cnsts.add(linExpr == 0);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    // Flow about delivery nodes; only when the warehouse visited
    //
    for (int k: prob->K_ae[a][e]) {
        linExpr.clear();
        sprintf(buf, "tFC(%d)(%d)(%d)", a, e, k);
        for (int j: prob->N_ae[a][e]) {
            linExpr += x_aeij[a][e][prob->n_k[k]][j];
        }
        for (int j: prob->N_ae[a][e]) {
            linExpr -= x_aeij[a][e][j][prob->h_k[k]];
        }
        cnsts.add(linExpr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    // Flow conservation
    //
    for (int i: prob->PD_ae[a][e]) {
        linExpr.clear();
        sprintf(buf, "FC_1(%d)(%d)(%d)", a, e, i);
        for (int j: prob->N_ae[a][e]) {
            linExpr += x_aeij[a][e][i][j];
        }
        cnsts.add(linExpr <= 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
        //
        sprintf(buf, "FC(%d)(%d)(%d)", a, e, i);
        for (int j: prob->N_ae[a][e]) {
            linExpr -= x_aeij[a][e][j][i];
        }
        cnsts.add(linExpr == 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    linExpr.end();
    baseModel->add(cnsts);
}

void Base::def_AT_cnsts_aeGiven(int a, int e) {
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
    cnsts.add(u_aei[a][e][prob->o_ae[a][e]] == prob->al_i[prob->o_ae[a][e]]);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    // Arrival time calculation
    //
    for (int i: prob->N_ae[a][e]) {
        for (int j: prob->N_ae[a][e]) {
            if (i == prob->d_ae[a][e] && j == prob->o_ae[a][e]) {
                continue;
            }
            linExpr.clear();
            sprintf(buf, "AT(%d)(%d)(%d)(%d)", a, e, i, j);
            linExpr += u_aei[a][e][i] + prob->t_ij[i][j];
            linExpr -= u_aei[a][e][j] + prob->M * (1 - x_aeij[a][e][i][j]);
            cnsts.add(linExpr <= 0);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
    }
    //
    // Time Window
    //
    for (int i: prob->N_ae[a][e]) {
        sprintf(buf, "TW1(%d)(%d)(%d)", a, e, i);
        cnsts.add(prob->al_i[i] <= u_aei[a][e][i]);
        cnsts[cnsts.getSize() - 1].setName(buf);
        //
        sprintf(buf, "TW2(%d)(%d)(%d)", a, e, i);
        cnsts.add(u_aei[a][e][i] <= prob->be_i[i]);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    // Routine route preservation
    //
    for (int i: prob->R_ae[a][e]) {
        for (int j: prob->R_ae[a][e]) {
            linExpr.clear();
            sprintf(buf, "RR_P(%d)(%d)(%d)(%d)", a, e, i, j);
            linExpr += prob->c_aeij[a][e][i][j] * u_aei[a][e][i];
            linExpr -= u_aei[a][e][j];
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
        linExpr += u_aei[a][e][prob->h_k[k]];
        linExpr -= u_aei[a][e][prob->n_k[k]];
        linExpr -= prob->M;
        for (int j: prob->N_ae[a][e]) {
            linExpr += prob->M * x_aeij[a][e][prob->n_k[k]][j];
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
            linExpr += prob->v_k[k] * x_aeij[a][e][j][prob->n_k[k]];
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
            linExpr += prob->w_k[k] * x_aeij[a][e][j][prob->n_k[k]];
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
            linExpr += prob->t_ij[i][j] * x_aeij[a][e][i][j];
        }
    }
    cnsts.add(linExpr <= prob->u_ae[a][e]);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.end();
    baseModel->add(cnsts);
}

void Base::def_COM_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
    IloExpr linExpr(env);
    //
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                linExpr.clear();
                sprintf(buf, "CC(%d)(%d)(%d)", a, e, k);
                linExpr += y_ak[a][k];
                for (int j: prob->N_ae[a][e]) {
                    linExpr -= x_aeij[a][e][j][prob->n_k[k]];
                }
                linExpr -= z_aek[a][e][k];
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