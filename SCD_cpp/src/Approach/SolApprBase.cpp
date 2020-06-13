//
//  SolApprBase.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 13/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/SolApprBase.hpp"


//std::mutex sab_mtx;

IloNumVar** new_inv_ak(Problem *prob, IloEnv &env, char vType) {
    IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
    //
    char buf[DEFAULT_BUFFER_SIZE];
    IloNumVar **y_ak = new IloNumVar*[prob->A.size()];
    for (int a: prob->A) {
        y_ak[a] = new IloNumVar[prob->K.size()];
        for (int k: prob->K) {
            sprintf(buf, "y(%d)(%d)", a, k);
            y_ak[a][k] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
        }
    }
    return y_ak;
}

IloNumVar*** new_inv_aek(Problem *prob, IloEnv &env, char vType) {
    IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
    //
    char buf[DEFAULT_BUFFER_SIZE];
    IloNumVar ***z_aek = new IloNumVar**[prob->A.size()];
    for (int a: prob->A) {
        z_aek[a] = new IloNumVar*[prob->E_a[a].size()];
        for (int e: prob->E_a[a]) {
            z_aek[a][e] = new IloNumVar[prob->K.size()];
            for (int k: prob->K) {
                sprintf(buf, "z(%d)(%d)(%d)", a, e, k);
                z_aek[a][e][k] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
            }
        }
    }
    return z_aek;
}

IloNumVar**** new_inv_aeij(Problem *prob, IloEnv &env, char vType) {
    IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
    //
    char buf[DEFAULT_BUFFER_SIZE];
    IloNumVar ****x_aeij = new IloNumVar***[prob->A.size()];
    for (int a: prob->A) {
        x_aeij[a] = new IloNumVar**[prob->E_a[a].size()];
        for (int e: prob->E_a[a]) {
            x_aeij[a][e] = new IloNumVar*[prob->cN.size()];
            for (int i = 0; i < prob->cN.size(); i++) {
                x_aeij[a][e][i] = new IloNumVar[prob->cN.size()];
            }
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    sprintf(buf, "x(%d)(%d)(%d)(%d)", a, e, i, j);
                    x_aeij[a][e][i][j] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
                }
            }
            
        }
    }
    return x_aeij;
}

IloNumVar*** new_inv_aei(Problem *prob, IloEnv &env) {
    char buf[DEFAULT_BUFFER_SIZE];
    IloNumVar ***u_aei = new IloNumVar**[prob->A.size()];
    for (int a: prob->A) {
        u_aei[a] = new IloNumVar*[prob->E_a[a].size()];
        for (int e: prob->E_a[a]) {
            u_aei[a][e] = new IloNumVar[prob->cN.size()];
            for (int i: prob->N_ae[a][e]) {
                sprintf(buf, "u(%d)(%d)(%d)", a, e, i);
                u_aei[a][e][i] = IloNumVar(env, 0.0, DBL_MAX, ILOFLOAT, buf);
            }
        }
    }
    return u_aei;
}

void delete_inv_ak(Problem *prob, IloNumVar **inv_ak) {
    for (int a: prob->A) {
        delete [] inv_ak[a];
    }
    delete [] inv_ak;
}

void delete_inv_aek(Problem *prob, IloNumVar ***inv_aek) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            delete [] inv_aek[a][e];
        }
        delete [] inv_aek[a];
    }
    delete [] inv_aek;
}

void delete_inv_aeij(Problem *prob, IloNumVar ****inv_aeij) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i = 0; i < prob->cN.size(); i++) {
                delete [] inv_aeij[a][e][i];
            }
            delete [] inv_aeij[a][e];
        }
        delete [] inv_aeij[a];
    }
    delete [] inv_aeij;
}

void delete_inv_aei(Problem *prob, IloNumVar ***inv_aei) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            delete [] inv_aei[a][e];
        }
        delete [] inv_aei[a];
    }
    delete [] inv_aei;
}

void BaseMM::build_baseModel() {
    def_preprocessing();
    //
    def_objF();
    def_ETA_cnsts(prob, env, y_ak, z_aek, baseModel);
    def_RUT_cnsts();
    def_COM_cnsts();
}
void BaseMM::def_preprocessing() {
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
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
        cnsts.add(x_aeij[a][e][i][i] == 0);
        sprintf(buf, "pP(%d)(%d)(%d)(%d)", a, e, i, j);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    baseModel->add(cnsts);
}

void BaseMM::def_objF() {
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

void BaseMM::def_RUT_cnsts() {
//    auto build_rutConsts = [](BaseMM *mm, int a, int e) {
//        mm->def_FC_cnsts_aeGiven(a, e);
//        mm->def_AT_cnsts_aeGiven(a, e);
//    };
    
    // Routing constraints
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            def_FC_cnsts_aeGiven(a, e);
            def_AT_cnsts_aeGiven(a, e);
        }
    }
    
//    std::vector<std::shared_future<void>> threadJobs;
//    for (int a : prob->A) {
//        for (int e: prob->E_a[a]) {
//            std::shared_future<void> returnValue = (std::shared_future<void>) routerPool.push(build_rutConsts, this, a, e);
//            threadJobs.push_back(returnValue);
//        }
//    }
//    std::for_each(threadJobs.begin(), threadJobs.end(), [](std::shared_future<void> x){ x.get(); });
    
}

void BaseMM::def_FC_cnsts_aeGiven(int a, int e) {
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
    for (int i: prob->N_ae[a][e]) {
        sprintf(buf, "NS(%d)(%d)(%d)", a, e, i);  // No Self Flow; tightening bounds
        cnsts.add(x_aeij[a][e][i][i] == 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
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
//    std::unique_lock<std::mutex> lock(sab_mtx);
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
//    std::unique_lock<std::mutex> lock(sab_mtx);
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
