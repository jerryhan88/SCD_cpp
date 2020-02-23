//
//  ILP.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 17/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "Base.hpp"

ILP::ILP(Problem *prob, std::string logPath, TimeTracker *tt, std::string lpPath) : Base(prob, logPath, tt) {
    this->lpPath = lpPath;
    cplexModel = new IloModel(env);
    //
    def_dvs();
    build_baseModel();
    //
    cplex = new IloCplex(*cplexModel);
}

Solution* ILP::solve() {
    cplex->exportModel(lpPath.c_str());
    
    
    Solution *sol = new Solution(prob);
    cplex->solve();
    if (cplex->getStatus() == IloAlgorithm::Infeasible) {
        cplex->exportModel(lpPath.c_str());
    }
    sol->objV = cplex->getObjValue();
    sol->gap = cplex->getMIPRelativeGap();
    sol->cpuT = tt->get_elipsedTimeCPU();
    sol->wallT = tt->get_elipsedTimeWall();
    char note[2048];
    sprintf(note,
            "\"{numRows: %ld,numCols: %ld}\"",
            cplex->getNrows(),
            cplex->getNcols());
    sol->note = std::string(note);
    //
    sol->alloMem4dvs();
    size_t numNodes = prob->locID.size();
    for (int a: prob->A) {
        for (int k: prob->K) {
            sol->y_ak[a][k] = cplex->getValue(y_ak[a][k]);
        }
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                sol->z_aek[a][e][k] = cplex->getValue(z_aek[a][e][k]);
            }
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    sol->x_aeij[a][e][i][j] = cplex->getValue(x_aeij[a][e][i][j]);
                }
                sol->u_aei[a][e][i] = cplex->getValue(u_aei[a][e][i]);
            }
        }
    }
    return sol;
}

void ILP::def_dvs() {
    size_t numAgents, numTasks, numNodes;
    numAgents = prob->A.size();
    numTasks = prob->K.size();
    numNodes = prob->locID.size();
    //
    y_ak = new IloNumVar*[numAgents];
    z_aek = new IloNumVar**[numAgents];
    x_aeij = new IloNumVar***[numAgents];
    u_aei = new IloNumVar**[numAgents];
    for (int a: prob->A) {
        y_ak[a] = new IloNumVar[numTasks];
        size_t numRR = prob->E_a[a].size();
        z_aek[a] = new IloNumVar*[numRR];
        x_aeij[a] = new IloNumVar**[numRR];
        u_aei[a] = new IloNumVar*[numRR];
        for (int e: prob->E_a[a]) {
            z_aek[a][e] = new IloNumVar[numTasks];
            x_aeij[a][e] = new IloNumVar*[numNodes];
            u_aei[a][e] = new IloNumVar[numNodes];
            for (int i = 0; i < numNodes; i++) {
                x_aeij[a][e][i] = new IloNumVar[numNodes];
            }
        }
    }
    char buf[DEFAULT_BUFFER_SIZE];
    for (int a: prob->A) {
        for (int k: prob->K) {
            sprintf(buf, "y(%d)(%d)", a, k);
            y_ak[a][k] = IloNumVar(env, 0.0, 1.0, ILOINT, buf);
        }
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                sprintf(buf, "z(%d)(%d)(%d)", a, e, k);
                z_aek[a][e][k] = IloNumVar(env, 0.0, 1.0, ILOINT, buf);
                
            }
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    sprintf(buf, "x(%d)(%d)(%d)(%d)", a, e, i, j);
                    x_aeij[a][e][i][j] = IloNumVar(env, 0.0, 1.0, ILOINT, buf);
                }
                sprintf(buf, "u(%d)(%d)(%d)", a, e, i);
                u_aei[a][e][i] = IloNumVar(env, 0.0, DBL_MAX, ILOFLOAT, buf);
            }
        }
    }
}

void ILP::build_baseModel() {
    def_objF();
    def_ETA_cnsts();
    def_RUT_cnsts();
    def_COM_cnsts();
}

void ILP::def_COM_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                linExpr.clear();
                linExpr += y_ak[a][k];
                for (int j: prob->N_ae[a][e]) {
                    linExpr -= x_aeij[a][e][j][prob->n_k[k]];
                }
                linExpr -= z_aek[a][e][k];
                cnsts.add(linExpr <= 0);
            }
        }
    }
    //
    cplexModel->add(cnsts);
}

void ILP::def_FC_cnsts_aeGiven(int a, int e) {
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    sprintf(buf, "s0_%d_%d", a, e);
    int ae_o = prob->locID[buf];
    sprintf(buf, "s%d_%d_%d", (int) prob->S_ae[a][e].size(), a, e);
    int ae_d = prob->locID[buf];
    //
    // Initiate flow
    //
    linExpr.clear();
    for (int j: prob->N_ae[a][e]) {
        linExpr += x_aeij[a][e][ae_o][j];
    }
    cnsts.add(linExpr <= 1);
    //
    linExpr.clear();
    for (int j: prob->N_ae[a][e]) {
        linExpr += x_aeij[a][e][j][ae_d];
    }
    cnsts.add(linExpr <= 1);
    //
    for (int i: prob->S_ae[a][e]) {
        if (i == ae_o || i == ae_d) continue;
        linExpr.clear();
        for (int j: prob->N_ae[a][e]) {
            if (j == i) continue;
            linExpr += x_aeij[a][e][i][j];
        }
        cnsts.add(linExpr <= 1);
        //
        linExpr.clear();
        for (int j: prob->N_ae[a][e]) {
            if (j == i) continue;
            linExpr += x_aeij[a][e][j][i];
        }
        cnsts.add(linExpr <= 1);
    }
    //
    // No flow
    //
    linExpr.clear();
    for (int j: prob->N_ae[a][e]) {
        linExpr += x_aeij[a][e][j][ae_o];
    }
    cnsts.add(linExpr <= 0);
    //
    linExpr.clear();
    for (int j: prob->N_ae[a][e]) {
        linExpr += x_aeij[a][e][ae_d][j];
    }
    cnsts.add(linExpr <= 0);
    //
    for (int k: prob->iF_ae[a][e]) {
        linExpr.clear();
        for (int j: prob->N_ae[a][e]) {
            linExpr += x_aeij[a][e][prob->n_k[k]][j];
        }
        cnsts.add(linExpr <= 0);
    }
    //
    linExpr.clear();
    for (int i: prob->N_ae[a][e]) {
        linExpr += x_aeij[a][e][i][i];
    }
    cnsts.add(linExpr <= 0);
    //
    // Flow about delivery nodes; only when the warehouse visited
    //
    for (int k: prob->F_ae[a][e]) {
        linExpr.clear();
        for (int j: prob->N_ae[a][e]) {
            linExpr += x_aeij[a][e][prob->n_k[k]][j];
        }
        for (int j: prob->N_ae[a][e]) {
            linExpr -= x_aeij[a][e][j][prob->h_k[k]];
        }
        cnsts.add(linExpr <= 0);
    }
    //
    // Flow conservation
    //
    for (int i: prob->N) {
        linExpr.clear();
        for (int j: prob->N_ae[a][e]) {
            linExpr += x_aeij[a][e][i][j];
        }
        for (int j: prob->N_ae[a][e]) {
            linExpr -= x_aeij[a][e][j][i];
        }
        cnsts.add(linExpr <= 0);
    }
    //
    cplexModel->add(cnsts);
}

void ILP::def_AT_cnsts_aeGiven(int a, int e) {
    //
    // Arrival time calculation
    //
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    // Time Window
    //
    
    for (int i: prob->N_ae[a][e]) {
        cnsts.add(u_aei[a][e][i] >= prob->al_i[i]);
        cnsts.add(u_aei[a][e][i] <= prob->be_i[i]);
    }
    //
    // Warehouse and Delivery Sequence
    //
    for (int k: prob->F_ae[a][e]) {
        linExpr.clear();
        linExpr = u_aei[a][e][prob->h_k[k]] - u_aei[a][e][prob->n_k[k]];
        cnsts.add(linExpr <= 0);
    }
    //
    // Routine route preservation
    //
    for (int i: prob->S_ae[a][e]) {
        for (int j: prob->S_ae[a][e]) {
            linExpr.clear();
            linExpr += prob->c_aeij[a][e][i][j] * u_aei[a][e][i];
            linExpr -= u_aei[a][e][j];
            cnsts.add(linExpr <= 0);
        }
    }
    //
    // Arrival time calculation
    //
    for (int i: prob->N_ae[a][e]) {
        for (int j: prob->N_ae[a][e]) {
            linExpr.clear();
            linExpr += u_aei[a][e][i];
            linExpr += prob->M * x_aeij[a][e][i][j];
            linExpr -= u_aei[a][e][j];
            cnsts.add(linExpr <= prob->M - prob->ga_i[i] - prob->t_ij[i][j]);
        }
    }
    //
    // Detour Limit
    //
    linExpr.clear();
    for (int i: prob->N_ae[a][e]) {
        for (int j: prob->N_ae[a][e]) {
            linExpr += (prob->ga_i[i] + prob->t_ij[i][j]) * x_aeij[a][e][i][j];
        }
    }
    cnsts.add(linExpr <= prob->l_ae[a][e] + prob->u_ae[a][e]);
    cplexModel->add(cnsts);
}

void ILP::def_RUT_cnsts() {
    // Routing constraints
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            def_FC_cnsts_aeGiven(a, e);
            def_AT_cnsts_aeGiven(a, e);
        }
    }
}

void ILP::def_ETA_cnsts() {
    // Evalutation of the Task Assignment
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    for (int k : prob->K) {
        linExpr.clear();
        for (int a : prob->A) {
            linExpr += y_ak[a][k];
        }
        cnsts.add(linExpr <= 1);
//        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    for (int a : prob->A) {
        linExpr.clear();
        for (int k : prob->K) {
            linExpr += prob->v_k[k] * y_ak[a][k];
        }
        cnsts.add(linExpr <= prob->v_a[a]);
    }
    //
    for (int a : prob->A) {
        linExpr.clear();
        for (int k : prob->K) {
            linExpr += prob->w_k[k] * y_ak[a][k];
        }
        cnsts.add(linExpr <= prob->w_a[a]);
    }
    //
    for (int k : prob->K) {
        for (int a : prob->A) {
            for (int e: prob->E_a[a]) {
                linExpr.clear();
                linExpr += z_aek[a][e][k];
                linExpr -= y_ak[a][k];
                cnsts.add(linExpr <= 0);
            }
        }
    }
    cplexModel->add(cnsts);
}

void ILP::def_objF() {
    IloExpr objF(env);
    for (int k: prob->K) {
        for (int a : prob->A) {
            objF += prob->r_k[k] * y_ak[a][k];
            for (int e: prob->E_a[a]) {
                objF -= prob->r_k[k] * prob->p_ae[a][e] * z_aek[a][e][k];
            }
        }
    }
    cplexModel->add(IloMaximize(env, objF));
    objF.end();
}
