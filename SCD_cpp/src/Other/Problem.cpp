//
//  Problem.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 16/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Problem.hpp"

#include <iostream>

double** new_dbl_ak(Problem *prob) {
    size_t numAgents, numTasks;
    numAgents = prob->A.size();
    numTasks = prob->K.size();
    //
    double **dbl_ak = new double*[numAgents];
    for (int a: prob->A) {
        dbl_ak[a] = new double[numTasks];
        for (int k: prob->K) {
            dbl_ak[a][k] = 0.0;
        }
    }
    return dbl_ak;
}

double*** new_dbl_aek(Problem *prob) {
    size_t numAgents, numTasks, numRR;
    numAgents = prob->A.size();
    numTasks = prob->K.size();
    //
    double ***dbl_aek = new double**[numAgents];
    for (int a: prob->A) {
        numRR = prob->E_a[a].size();
        dbl_aek[a] = new double*[numRR];
        for (int e: prob->E_a[a]) {
            dbl_aek[a][e] = new double[numTasks];
            for (int k: prob->K) {
                dbl_aek[a][e][k] = 0.0;
            }
        }
    }
    return dbl_aek;
}

double**** new_dbl_aeij(Problem *prob) {
    size_t numAgents, numTasks, numNodes, numRR;
    numAgents = prob->A.size();
    numTasks = prob->K.size();
    numNodes = prob->cN.size();
    //
    double ****dbl_aeij = new double***[numAgents];
    for (int a: prob->A) {
        numRR = prob->E_a[a].size();
        dbl_aeij[a] = new double**[numRR];
        for (int e: prob->E_a[a]) {
            dbl_aeij[a][e] = new double*[numNodes];
            for (int i = 0; i < numNodes; i++) {
                dbl_aeij[a][e][i] = new double[numNodes];
                for (int j = 0; j < numNodes; j++) {
                    dbl_aeij[a][e][i][j] = 0.0;
                }
            }
        }
    }
    return dbl_aeij;
}


double*** new_dbl_aei(Problem *prob) {
    size_t numAgents, numTasks, numNodes, numRR;
    numAgents = prob->A.size();
    numTasks = prob->K.size();
    numNodes = prob->cN.size();
    //
    double ***dbl_aei = new double**[numAgents];
    for (int a: prob->A) {
        numRR = prob->E_a[a].size();
        dbl_aei[a] = new double*[numRR];
        for (int e: prob->E_a[a]) {
            dbl_aei[a][e] = new double[numNodes];
            for (int i = 0; i < numNodes; i++) {
                dbl_aei[a][e][i] = 0.0;
            }
        }
    }
    return dbl_aei;
}

void delete_dbl_ak(Problem *prob, double **dbl_ak) {
    for (int a: prob->A) {
        delete [] dbl_ak[a];
    }
    delete [] dbl_ak;
}

void delete_dbl_aek(Problem *prob, double ***dbl_aek) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            delete [] dbl_aek[a][e];
        }
        delete [] dbl_aek[a];
    }
    delete [] dbl_aek;
}

void delete_dbl_aeij(Problem *prob, double ****dbl_aeij) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int i = 0; i < prob->cN.size(); i++) {
                delete [] dbl_aeij[a][e][i];
            }
            delete [] dbl_aeij[a][e];
        }
        delete [] dbl_aeij[a];
    }
    delete [] dbl_aeij;
}

void delete_dbl_aei(Problem *prob, double ***dbl_aei) {
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            delete [] dbl_aei[a][e];
        }
        delete [] dbl_aei[a];
    }
    delete [] dbl_aei;
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

