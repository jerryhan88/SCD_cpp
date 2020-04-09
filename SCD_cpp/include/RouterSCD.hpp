//
//  RouterSCD.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 25/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef RouterSCD_h
#define RouterSCD_h

#include <string>

#include <ilcplex/ilocplex.h>

#include "Problem.hpp"
//#include "BaseMM.hpp"

#include "ck_route/Router.hpp"     // from BnC_CPLEX
#include "ck_route/Other.hpp"     // from BnC_CPLEX
#include "ck_util/util.hpp"         // from util

void update_objF(rmm::RouteMM *rutMM, Problem *prob, int a, int e, double ***lrh_l_aek);
void run_rutAlgo(rmm::RouteMM *rutMM, double *rut_objV, double **rut_x_ij, double *rut_u_i);

class RouterSCD {
public:
    Problem *prob;
    TimeTracker *tt;
    rut::Problem *rutProb;
    int a, e;
    double ***lrh_l_aek;
    double **rut_x_ij, *rut_u_i;
    double rut_objV;
    //
    IloCplex *rutCplex;
    //
    RouterSCD(Problem *prob, TimeTracker *tt, int a, int e, double ***lrh_l_aek) {
        this->prob = prob;
        this->tt = tt;
        rutProb = prob->RP_ae[a][e];
        this->a = a;
        this->e = e;
        this->lrh_l_aek = lrh_l_aek;
        //
        rut_x_ij = new double*[rutProb->N.size()];
        rut_u_i = new double[rutProb->N.size()];
        for (int i: rutProb->N) {
            rut_x_ij[i] = new double[rutProb->N.size()];
            for (int j: rutProb->N) {
                rut_x_ij[i][j] = 0.0;
            }
            rut_u_i[i] = 0.0;
        }
    }
    ~RouterSCD() {
        for (int i: rutProb->N) {
            delete [] rut_x_ij[i];
        }
        delete [] rut_x_ij;
        delete [] rut_u_i;
    }
    virtual void update() {
        throw "Should override update()";
    }
    virtual void run() {
        throw "Should override run()";
    }
};

class RouterILP: public RouterSCD {
public:
    std::string rutLogPath;
    rmm::ILP *rutAlgo;
    std::vector<std::string> tempLogMsg;
    //
    RouterILP(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec,
              int a, int e,
              double ***lrh_l_aek,
              std::string rutLogPath, std::string rutLpLogPath): RouterSCD(prob, tt, a, e, lrh_l_aek) {
        this->rutLogPath = rutLogPath;
        rutAlgo = new rmm::ILP(rutProb, tt, time_limit_sec, rutLogPath, rutLpLogPath, true);
        rutAlgo->cplex->setWarning(rutAlgo->env.getNullStream());
        this->rutCplex = rutAlgo->cplex;
    }
    ~RouterILP() {
        delete rutAlgo;
    }
    //
    void update();
    void run();
};

class RouterBnC: public RouterSCD {
public:
    const std::vector<std::string> cut_names {"hSE", "hCA", "hRS", "hIP"};
    const IloCplex::CutManagement cutManagerType = IloCplex::UseCutForce;
    const IloBool isLocalCutAdd = IloFalse;
    //
    rmm::BnC *rutAlgo;
    bool turOnCutPool;
    //
    RouterBnC(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec,
              int a, int e,
              double ***lrh_l_aek,
              std::string rutLogPath, std::string rutLpLogPath,
              bool turOnCutPool): RouterSCD(prob, tt, a, e, lrh_l_aek) {
        rutAlgo = new rmm::BnC(rutProb, tt, time_limit_sec, rutLogPath, rutLpLogPath, cut_names, true, cutManagerType, isLocalCutAdd);
        rutAlgo->cplex->setWarning(rutAlgo->env.getNullStream());
        this->turOnCutPool = turOnCutPool;
        this->rutCplex = rutAlgo->cplex;
    }
    //
    void update();
    void run();
};

class RouterGH: public RouterSCD {
public:
    rgh::InsertionHeuristic *rutAlgo;
    //
    RouterGH(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec,
            int a, int e,
            double ***lrh_l_aek,
            std::string rutLogPath): RouterSCD(prob, tt, a, e, lrh_l_aek) {
        rutAlgo = new rgh::InsertionHeuristic(rutProb, tt, time_limit_sec, rutLogPath);
    }
    //
    void update();
    void run();
};


#endif /* RouterSCD_h */
