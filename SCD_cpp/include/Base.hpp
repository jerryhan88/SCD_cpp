//
//  Base.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 17/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef Base_hpp
#define Base_hpp

#include <cfloat>
#include <math.h>

#include <ilcplex/ilocplex.h>

#include "Problem.hpp"
#include "Solution.hpp"

#include "ck_route/CutBase.hpp"     // from BnC_CPLEX
#include "ck_route/RouteMM.hpp"     // from BnC_CPLEX
#include "ck_util/util.hpp"         // from util


#define DEFAULT_BUFFER_SIZE 2048

class Base {
public:
    Problem *prob;
    std::string logPath, lpPath;
    TimeTracker *tt;
    //
    IloEnv env;
    IloModel *baseModel;
    IloCplex *baseCplex;
    IloNumVar **y_ak, ***z_aek, ****x_aeij, ***u_aei;
    IloRangeArray *COM_cnsts;
    long ***COM_cnsts_index;
    //
    Base(Problem *prob, std::string logPath, std::string lpPath, TimeTracker *tt, char vType) {
        this->prob = prob;
        this->logPath = logPath;
        this->lpPath = lpPath;
        this->tt = tt;
        //
        y_ak = Base::gen_y_ak(prob, env, vType);
        z_aek = Base::gen_z_aek(prob, env, vType);
        x_aeij = Base::gen_x_aeij(prob, env, vType);
        u_aei = Base::gen_u_aei(prob, env);
        COM_cnsts = new IloRangeArray(env);
        COM_cnsts_index = new long**[prob->A.size()];
        for (int a : prob->A) {
            COM_cnsts_index[a] = new long*[prob->E_a[a].size()];
            for (int e: prob->E_a[a]) {
                COM_cnsts_index[a][e] = new long[prob->K.size()];
            }
        }
        //
        baseModel = new IloModel(env);
        build_baseModel();
        baseCplex = new IloCplex(*baseModel);
        baseCplex->setOut(env.getNullStream());
    }
    ~Base() {
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int i = 0; i < prob->cN.size(); i++) {
                    delete [] x_aeij[a][e][i];
                }
                delete [] z_aek[a][e]; delete [] x_aeij[a][e]; delete [] u_aei[a][e];
            }
            delete [] y_ak[a]; delete [] z_aek[a]; delete [] x_aeij[a]; delete [] u_aei[a];
        }
        delete [] y_ak; delete [] z_aek; delete [] x_aeij; delete [] u_aei;
        //
        delete baseCplex;
        delete baseModel;
    }
    virtual Solution* solve() {
        throw "Should override solve()";
    }
    static IloNumVar** gen_y_ak(Problem *prob, IloEnv &env, char vType) {
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
    
    static IloNumVar*** gen_z_aek(Problem *prob, IloEnv &env, char vType) {
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
    
    static IloNumVar**** gen_x_aeij(Problem *prob, IloEnv &env, char vType) {
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
    
    static IloNumVar*** gen_u_aei(Problem *prob, IloEnv &env) {
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
    
    static void def_ETA_cnsts(Problem *prob, IloEnv &env, IloNumVar **y_ak, IloNumVar ***z_aek, IloModel *model) {
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
                linExpr += y_ak[a][k];
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
                        sprintf(buf, "IA(%d)(%d)(%d)", a, e, k);  // Infeasible Assignment
                        linExpr += y_ak[a][k];
                        cnsts.add(linExpr == 0);
                        cnsts[cnsts.getSize() - 1].setName(buf);
                    }
                }
            }
        }
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int k: prob->K) {
                    linExpr.clear();
                    sprintf(buf, "TAC(%d)(%d)(%d)", a, e, k);  // Task Accomplishment
                    linExpr += z_aek[a][e][k];
                    linExpr -= y_ak[a][k];
                    cnsts.add(linExpr <= 0);
                    cnsts[cnsts.getSize() - 1].setName(buf);
                }
            }
        }
        model->add(cnsts);
    }
    
protected:
    void build_baseModel();
    void def_objF();
    void def_RUT_cnsts();
    void def_FC_cnsts_aeGiven(int a, int e);
    void def_AT_cnsts_aeGiven(int a, int e);
    void def_COM_cnsts();
};
    
class ILP : public Base {
public:
    ILP(Problem *prob, std::string logPath, std::string lpPath, TimeTracker *tt) : Base(prob, logPath, lpPath, tt, 'I') {
        this->lpPath = lpPath;
    }
    ~ILP() {
    }
    Solution* solve();
};
    
class LRH : public Base {
public:
    std::string _ROUTER;
    double STEP_DECREASE_RATE, DUAL_GAP_LIMIT;
    unsigned int NO_IMPROVEMENT_LIMIT, NUM_ITER_LIMIT;
    unsigned long TIME_LIMIT_SEC;
    bool turnOnInitSol, turOnCutPool;
    //
    double L1_V, L2_V, L_V, L_V_star, F_V, F_V_star;
    double ut;
    int numIters, noUpdateCounter_L;
    double **lrh_y_ak, ***lrh_z_aek, ****lrh_x_aeij, ***lrh_u_aei, ***lrh_l_aek;
    //
    IloModel *etaModel; // Evaluation on Task Assingment
    IloModel *pexModel;  // Primal Extraction
    IloCplex *etaCplex, *pexCplex;
    IloNumVar **eta_y_ak, ***eta_z_aek;
    IloNumVar **pex_y_ak, ***pex_z_aek;
    IloRangeArray *pex_COM_cnsts;
    long ***pex_COM_cnsts_index;
    std::vector<std::vector<rmm::RouteMM*>> rutMMs;
    std::vector<std::vector<rmm::BnC*>> rutBnCs;
    std::vector<std::vector<double**>> rut_x_ij;
    std::vector<std::vector<double*>> rut_u_i;
    //
    Solution *bestSol;
    //
    LRH(Problem *prob, std::string logPath, std::string lpPath, TimeTracker *tt,
        std::string _router,
        double dual_gap_limit, unsigned int num_iter_limit, unsigned int no_improvement_limit,
        unsigned long time_limit_sec,
        bool turnOnInitSol, bool turOnCutPool) : Base(prob, logPath, lpPath, tt, 'C') {
        _ROUTER = _router;
        DUAL_GAP_LIMIT = dual_gap_limit;
        NUM_ITER_LIMIT = num_iter_limit;
        NO_IMPROVEMENT_LIMIT = no_improvement_limit;
        TIME_LIMIT_SEC = time_limit_sec;
        this->turnOnInitSol = turnOnInitSol;
        this->turOnCutPool = turOnCutPool;
        L_V_star = DBL_MAX; F_V_star = -DBL_MAX;
        STEP_DECREASE_RATE = 0.5;
        ut = 2.0;
        noUpdateCounter_L = 0;
        numIters = 0;
        //
        lrh_y_ak = new double*[prob->A.size()];
        lrh_z_aek = new double**[prob->A.size()];
        lrh_x_aeij = new double***[prob->A.size()];
        lrh_u_aei = new double**[prob->A.size()];
        lrh_l_aek = new double**[prob->A.size()];
        for (int a: prob->A) {
            lrh_y_ak[a] = new double[prob->K.size()];
            for (int k: prob->K) {
                lrh_y_ak[a][k] = 0.0;
            }
            lrh_z_aek[a] = new double*[prob->E_a[a].size()];
            lrh_x_aeij[a] = new double**[prob->E_a[a].size()];
            lrh_u_aei[a] = new double*[prob->E_a[a].size()];
            lrh_l_aek[a] = new double*[prob->E_a[a].size()];
            for (int e: prob->E_a[a]) {
                lrh_z_aek[a][e] = new double[prob->K.size()];
                lrh_x_aeij[a][e] = new double*[prob->cN.size()];
                lrh_u_aei[a][e] = new double[prob->cN.size()];
                lrh_l_aek[a][e] = new double[prob->K.size()];
                for (int k: prob->K) {
                    lrh_z_aek[a][e][k] = 0.0;
                    lrh_l_aek[a][e][k] = 0.0;
                }
                for (int i = 0; i < prob->cN.size(); i++) {
                    lrh_x_aeij[a][e][i] = new double[prob->cN.size()];
                    lrh_u_aei[a][e][i] = 0.0;
                    for (int j = 0; j < prob->cN.size(); j++) {
                        lrh_x_aeij[a][e][i][j] = 0.0;
                    }
                }
            }
        }
        pex_COM_cnsts = new IloRangeArray(env);
        pex_COM_cnsts_index = new long**[prob->A.size()];
        for (int a : prob->A) {
            pex_COM_cnsts_index[a] = new long*[prob->E_a[a].size()];
            for (int e: prob->E_a[a]) {
                pex_COM_cnsts_index[a][e] = new long[prob->K.size()];
            }
        }
        //
        build_etaModel();
        build_pexModel();
        build_rutModels();
        //
        bestSol = nullptr;
        //
        if (logPath != "") {
            std::string _header("wallT,cpuT,Iteration,Function,L_V*,L_V,F_V*,F_V,Note");
            char header[_header.size() + 1];
            std::strcpy(header, _header.c_str());
            createCSV(logPath, header);
        }
    }
    ~LRH() {
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int i = 0; i < prob->cN.size(); i++) {
                    delete [] lrh_x_aeij[a][e][i];
                }
                delete [] lrh_z_aek[a][e]; delete [] lrh_x_aeij[a][e]; delete [] lrh_l_aek[a][e];
                delete [] eta_z_aek[a][e];
                delete [] pex_z_aek[a][e];
                delete [] pex_COM_cnsts_index[a][e];
                delete rutMMs[a][e];
                if (turnOnInitSol) {
                    for (int i: prob->RP_ae[a][e]->N) {
                        delete [] rut_x_ij[a][e][i];
                    }
                    delete [] rut_x_ij[a][e];
                    delete [] rut_u_i[a][e];
                }
            }
            delete [] lrh_y_ak[a]; delete [] lrh_z_aek[a]; delete [] lrh_x_aeij[a]; delete [] lrh_l_aek[a];
            delete [] eta_y_ak[a]; delete [] eta_z_aek[a];
            delete [] pex_y_ak[a]; delete [] pex_z_aek[a];
            delete [] pex_COM_cnsts_index[a];
        }
        delete [] lrh_y_ak; delete [] lrh_z_aek; delete [] lrh_x_aeij; delete [] lrh_l_aek;
        delete [] eta_y_ak; delete [] eta_z_aek;
        delete [] pex_y_ak; delete [] pex_z_aek;
        delete [] pex_COM_cnsts_index;
        //
        delete etaCplex; delete etaModel;
        delete pexCplex; delete pexModel;
        delete pex_COM_cnsts;
        env.end();
        rutMMs.clear(); rutBnCs.clear(); rut_x_ij.clear(); rut_u_i.clear();
    }
    //
    Solution* solve();
private:
    void solve_dualProblem();
    void updateLMs();
    //
    void build_etaModel();
    void update_etaModel();
    void solve_etaModel();
    //
    void build_rutModels();
    void update_rutMM(int a, int e);
    void solve_rutModels();
    //
    void build_pexModel();
    void update_pexModel();
    void solve_pexModel();
    //
    void logging(std::string indicator,
                 std::string _L_V_star, std::string _L_V, std::string _F_V_star, std::string _F_V,
                 std::string note);

    
};
#endif /* Base_hpp */
