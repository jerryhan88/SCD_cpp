//
//  Problem.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 16/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Problem.hpp"


IloNumVar** new_ak_inv(Problem *prob, IloEnv &env, char vType) {
    IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
    //
    char buf[2048];
    IloNumVar **ak_ary = new IloNumVar*[prob->A.size()];
    for (int a: prob->A) {
        ak_ary[a] = new IloNumVar[prob->K.size()];
        for (int k: prob->K) {
            sprintf(buf, "y(%d)(%d)", a, k);
            ak_ary[a][k] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
        }
    }
    return ak_ary;
}
void new_ak_inv(Problem *prob, IloEnv &env, char vType, std::map<iiTup, IloNumVar> &ak_map) {
    IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
    //
    char buf[2048];
    for (int a: prob->A) {
        for (int k: prob->K) {
            sprintf(buf, "y(%d)(%d)", a, k);
            ak_map[std::make_tuple(a, k)] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
        }
    }
}
void del_ak_inv(Problem *prob, IloNumVar **ak_ary) {
    for (int a: prob->A) {
        delete [] ak_ary[a];
    }
    delete [] ak_ary;
}

IloNumVar*** new_aek_inv(Problem *prob, IloEnv &env, char vType) {
    IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
    //
    char buf[2048];
    IloNumVar ***aek_ary = new IloNumVar**[prob->A.size()];
    for (int a: prob->A) {
        aek_ary[a] = new IloNumVar*[prob->E_a[a].size()];
        for (int e: prob->E_a[a]) {
            aek_ary[a][e] = new IloNumVar[prob->K.size()];
            for (int k: prob->K) {
                sprintf(buf, "z(%d)(%d)(%d)", a, e, k);
                aek_ary[a][e][k] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
            }
        }
    }
    return aek_ary;
}
void new_aek_inv(Problem *prob, IloEnv &env, char vType, std::map<iiiTup, IloNumVar> &aek_map) {
    IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
    //
    char buf[2048];
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                sprintf(buf, "z(%d)(%d)(%d)", a, e, k);
                aek_map[std::make_tuple(a, e, k)] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
            }
        }
    }
}
void del_aek_inv(Problem *prob, IloNumVar ***aek_ary) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            delete [] aek_ary[a][e];
        }
        delete [] aek_ary[a];
    }
    delete [] aek_ary;
}

IloNumVar**** new_aeij_inv(Problem *prob, IloEnv &env, char vType, bool ****bool_x_aeij) {
    IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
    //
    char buf[2048];
    IloNumVar ****aeij_ary = new IloNumVar***[prob->A.size()];
    for (int a: prob->A) {
        aeij_ary[a] = new IloNumVar**[prob->E_a[a].size()];
        for (int e: prob->E_a[a]) {
            aeij_ary[a][e] = new IloNumVar*[prob->cN.size()];
            for (int i = 0; i < prob->cN.size(); i++) {
                aeij_ary[a][e][i] = new IloNumVar[prob->cN.size()];
            }
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    if (!bool_x_aeij[a][e][i][j]) {
                        continue;
                    }
                    sprintf(buf, "x(%d)(%d)(%d)(%d)", a, e, i, j);
                    aeij_ary[a][e][i][j] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
                }
            }
            
        }
    }
    return aeij_ary;
}
void new_aeij_inv(Problem *prob, IloEnv &env, char vType, bool ****bool_x_aeij,
                  std::map<iiiiTup, IloNumVar> &aeij_map) {
    IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
    //
    char buf[2048];
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    if (!bool_x_aeij[a][e][i][j]) {
                        continue;
                    }
                    sprintf(buf, "x(%d)(%d)(%d)(%d)", a, e, i, j);
                    aeij_map[std::make_tuple(a, e, i, j)] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
                }
            }
        }
    }
}
void del_aeij_inv(Problem *prob, IloNumVar ****aeij_ary) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i = 0; i < prob->cN.size(); i++) {
                delete [] aeij_ary[a][e][i];
            }
            delete [] aeij_ary[a][e];
        }
        delete [] aeij_ary[a];
    }
    delete [] aeij_ary;
}

IloNumVar*** new_aei_inv(Problem *prob, IloEnv &env) {
    char buf[2048];
    IloNumVar ***aei_ary = new IloNumVar**[prob->A.size()];
    for (int a: prob->A) {
        aei_ary[a] = new IloNumVar*[prob->E_a[a].size()];
        for (int e: prob->E_a[a]) {
            aei_ary[a][e] = new IloNumVar[prob->cN.size()];
            for (int i: prob->N_ae[a][e]) {
                sprintf(buf, "u(%d)(%d)(%d)", a, e, i);
                aei_ary[a][e][i] = IloNumVar(env, 0.0, DBL_MAX, ILOFLOAT, buf);
            }
        }
    }
    return aei_ary;
}
void new_aei_inv(Problem *prob, IloEnv &env, std::map<iiiTup, IloNumVar> &aei_map) {
    char buf[2048];
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i: prob->N_ae[a][e]) {
                sprintf(buf, "u(%d)(%d)(%d)", a, e, i);
                aei_map[std::make_tuple(a, e, i)] = IloNumVar(env, 0.0, DBL_MAX, ILOFLOAT, buf);
            }
        }
    }
}
void del_aei_inv(Problem *prob, IloNumVar ***aei_ary) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            delete [] aei_ary[a][e];
        }
        delete [] aei_ary[a];
    }
    delete [] aei_ary;
}


double** new_ak_dbl(Problem *prob) {
    size_t numAgents, numTasks;
    numAgents = prob->A.size();
    numTasks = prob->K.size();
    //
    double **ak_ary = new double*[numAgents];
    for (int a: prob->A) {
        ak_ary[a] = new double[numTasks];
        for (int k: prob->K) {
            ak_ary[a][k] = 0.0;
        }
    }
    return ak_ary;
}
void new_ak_dbl(Problem *prob, std::map<iiTup, double> &ak_map) {
    for (int a: prob->A) {
        for (int k: prob->K) {
            ak_map[std::make_tuple(a, k)] = 0.0;
        }
    }
}
void del_ak_dbl(Problem *prob, double **ak_ary) {
    for (int a: prob->A) {
        delete [] ak_ary[a];
    }
    delete [] ak_ary;
}

double*** new_aek_dbl(Problem *prob) {
    size_t numAgents, numTasks, numRR;
    numAgents = prob->A.size();
    numTasks = prob->K.size();
    //
    double ***aek_ary = new double**[numAgents];
    for (int a: prob->A) {
        numRR = prob->E_a[a].size();
        aek_ary[a] = new double*[numRR];
        for (int e: prob->E_a[a]) {
            aek_ary[a][e] = new double[numTasks];
            for (int k: prob->K) {
                aek_ary[a][e][k] = 0.0;
            }
        }
    }
    return aek_ary;
}
void new_aek_dbl(Problem *prob, std::map<iiiTup, double> &aek_map) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                aek_map[std::make_tuple(a, e, k)] = 0.0;
            }
        }
    }
}
void del_aek_dbl(Problem *prob, double ***aek_ary) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            delete [] aek_ary[a][e];
        }
        delete [] aek_ary[a];
    }
    delete [] aek_ary;
}

double**** new_aeij_dbl(Problem *prob) {
    size_t numAgents, numNodes, numRR;
    numAgents = prob->A.size();
    numNodes = prob->cN.size();
    //
    double ****aeij_ary = new double***[numAgents];
    for (int a: prob->A) {
        numRR = prob->E_a[a].size();
        aeij_ary[a] = new double**[numRR];
        for (int e: prob->E_a[a]) {
            aeij_ary[a][e] = new double*[numNodes];
            for (int i = 0; i < numNodes; i++) {
                aeij_ary[a][e][i] = new double[numNodes];
                for (int j = 0; j < numNodes; j++) {
                    aeij_ary[a][e][i][j] = 0.0;
                }
            }
        }
    }
    return aeij_ary;
}
void new_aeij_dbl(Problem *prob, bool ****bool_x_aeij, std::map<iiiiTup, double> &aeij_map) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    if (!bool_x_aeij[a][e][i][j]) {
                        continue;
                    }
                    aeij_map[std::make_tuple(a, e, i, j)] = 0.0;
                }
            }
        }
    }
}
void del_aeij_dbl(Problem *prob, double ****aeij_ary) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i = 0; i < prob->cN.size(); i++) {
                delete [] aeij_ary[a][e][i];
            }
            delete [] aeij_ary[a][e];
        }
        delete [] aeij_ary[a];
    }
    delete [] aeij_ary;
}

double*** new_aei_dbl(Problem *prob) {
    size_t numAgents, numNodes, numRR;
    numAgents = prob->A.size();
    numNodes = prob->cN.size();
    //
    double ***aei_ary = new double**[numAgents];
    for (int a: prob->A) {
        numRR = prob->E_a[a].size();
        aei_ary[a] = new double*[numRR];
        for (int e: prob->E_a[a]) {
            aei_ary[a][e] = new double[numNodes];
            for (int i = 0; i < numNodes; i++) {
                aei_ary[a][e][i] = 0.0;
            }
        }
    }
    return aei_ary;
}
void new_aei_dbl(Problem *prob, std::map<iiiTup, double> &aei_map) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i: prob->N_ae[a][e]) {
                aei_map[std::make_tuple(a, e, i)] = 0.0;
            }
        }
    }
}
void del_aei_dbl(Problem *prob, double ***aei_ary) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            delete [] aei_ary[a][e];
        }
        delete [] aei_ary[a];
    }
    delete [] aei_ary;
}

void Problem::gen_aeProbs() {
    //
    // Generate problem instances associated with routing problems by rearranging indices of parameters
    //
    for (int a: A) {
        std::vector<rut::Problem*> aRP;
        std::vector<std::map<int, int>> a_scdID_rutID, a_rutID_scdID;
        for (int e: E_a[a]) {
            //
            // Rearrange indices
            //
            std::map<int, int> scdID_rutID, rutID_scdID;
            int o = 0;
            int d = (int) P_ae[a][e].size() + (int) D_ae[a][e].size() + 1;
            rutID_scdID[o] = o_ae[a][e];
            scdID_rutID[o_ae[a][e]] = o;
            for (int scdID: P_ae[a][e]) {
                int rudID = (int) rutID_scdID.size();
                rutID_scdID[rudID] = scdID;
                scdID_rutID[scdID] = rudID;
            }
            for (int scdID: D_ae[a][e]) {
                int rudID = (int) rutID_scdID.size();
                rutID_scdID[rudID] = scdID;
                scdID_rutID[scdID] = rudID;
            }
            rutID_scdID[d] = d_ae[a][e];
            scdID_rutID[d_ae[a][e]] = d;
            for (int scdID: R_ae[a][e]) {
                if (scdID == o_ae[a][e] || scdID == d_ae[a][e]) {
                    continue;
                }
                int rudID = (int) rutID_scdID.size();
                rutID_scdID[rudID] = scdID;
                scdID_rutID[scdID] = rudID;
            }
            a_scdID_rutID.push_back(scdID_rutID);
            a_rutID_scdID.push_back(rutID_scdID);
            //
            // Generate a problem instance associated with a routine route and its feasible tasks
            //
            rut::Problem *aeProb = new rut::Problem();
            std::vector<int> _aeK(K_ae[a][e].begin(), K_ae[a][e].end());
            for (int k = 0; k < _aeK.size(); k++) {
                aeProb->K.push_back(k);
            }
            for (int scdID: PD_ae[a][e]) {
                aeProb->PD.push_back(scdID_rutID[scdID]);
            }
            for (int scdID: N_ae[a][e]) {
                aeProb->N.push_back(scdID_rutID[scdID]);
            }
            for (int scdID: R_ae[a][e]) {
                aeProb->R.insert(scdID_rutID[scdID]);
                aeProb->_R.push_back(scdID_rutID[scdID]);
            }
            for (int i = 0; i < aeProb->_R.size() - 1; i++) {
                std::set<int> lS;
                lS.insert(aeProb->_R.begin() + i + 2, aeProb->_R.end());
                aeProb->LR.push_back(lS);
            }
            for (int scdID: P_ae[a][e]) {
                aeProb->P.insert(scdID_rutID[scdID]);
            }
            std::vector<int> _D;
            for (int scdID: D_ae[a][e]) {
                aeProb->D.insert(scdID_rutID[scdID]);
                _D.push_back(scdID_rutID[scdID]);
            }
            aeProb->o = o;
            aeProb->d = d;
            aeProb->r_k = new double[aeProb->K.size()];
            aeProb->v_k = new double[aeProb->K.size()];
            aeProb->w_k = new double[aeProb->K.size()];
            aeProb->h_k = new int[aeProb->K.size()];
            aeProb->n_k = new int[aeProb->K.size()];
            for (int k: aeProb->K) {
                aeProb->r_k[k] = r_k[_aeK[k]];
                aeProb->v_k[k] = v_k[_aeK[k]];
                aeProb->w_k[k] = w_k[_aeK[k]];
                aeProb->h_k[k] = scdID_rutID[h_k[_aeK[k]]];
                aeProb->n_k[k] = scdID_rutID[n_k[_aeK[k]]];
            }
            //
            size_t numNodes = scdID_rutID.size();
            aeProb->t_ij = new double*[numNodes];
            aeProb->c_ij = new int*[numNodes];
            for (int i = 0; i < numNodes; i++) {
                aeProb->t_ij[i] = new double[numNodes];
                aeProb->c_ij[i] = new int[numNodes];
            }
            aeProb->al_i = new double[numNodes];
            aeProb->be_i = new double[numNodes];
            for (int rutID0: aeProb->N) {
                int scdID0 = rutID_scdID[rutID0];
                aeProb->al_i[rutID0] = al_i[scdID0];
                aeProb->be_i[rutID0] = be_i[scdID0];
                assert(0 <= aeProb->be_i[rutID0]);
                for (int rutID1: aeProb->N) {
                    int scdID1 = rutID_scdID[rutID1];
                    aeProb->t_ij[rutID0][rutID1] = t_ij[scdID0][scdID1];
                    aeProb->c_ij[rutID0][rutID1] = c_aeij[a][e][scdID0][scdID1];
                }
            }
            aeProb->M = M;
            aeProb->bv = v_a[a];
            aeProb->bw = w_a[a];
            aeProb->bu = u_ae[a][e];
            //
            aeProb->v_i = new double[numNodes];
            aeProb->w_i = new double[numNodes];
            for (int n0: aeProb->R) {
                aeProb->v_i[n0] = 0.0;
                aeProb->w_i[n0] = 0.0;
            }
            for (int n0: aeProb->P) {
                aeProb->v_i[n0] = 0.0;
                aeProb->w_i[n0] = 0.0;
            }
            for (int k0 = 0; k0 < _D.size(); k0++) {
                int rutID = _D[k0];
                int k1 = -1;
                for (int i = 0; i < _D.size(); i++) {
                    if (aeProb->n_k[i] == rutID) {
                        k1 = i;
                        break;
                    }
                }
                assert(k1 != -1);
                aeProb->v_i[rutID] = aeProb->v_k[k1];
                aeProb->w_i[rutID] = aeProb->w_k[k1];
            }
            aRP.push_back(aeProb);
        }
        RP_ae.push_back(aRP);
        scdID_rutID_ae.push_back(a_scdID_rutID);
        rutID_scdID_ae.push_back(a_rutID_scdID);
    }
}

