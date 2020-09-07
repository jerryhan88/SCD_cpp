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

#include "ck_route/Router.hpp"     // from BnC_CPLEX
#include "ck_route/Other.hpp"     // from BnC_CPLEX
#include "ck_util/util.hpp"         // from util

void update_objF(rmm::RouteMM *rutMM, Problem *prob, int a, int e, double ***lrh_l_ary);
void run_rutAlgo(rmm::RouteMM *rutMM, double *rut_objV, double **rut_x_ij, double *rut_u_i);

class RouterSCD {
public:
    Problem *prob;
    TimeTracker *tt;
    unsigned long time_limit_sec;
    rut::Problem *rutProb;
    int a, e;
    double ***lrh_l_ary;
    std::map<iiiTup, double> *lrh_l_map;
    double **rut_x_ij, *rut_u_i;
    double rut_objV;
    //
    IloCplex *rutCplex;
    //
    RouterSCD(Problem *prob, TimeTracker *tt,
              unsigned long time_limit_sec,
              int a, int e,
              double ***lrh_l_ary, std::map<iiiTup, double> *lrh_l_map) {
        this->prob = prob;
        this->tt = tt;
        this->time_limit_sec = time_limit_sec;
        rutProb = prob->RP_ae[a][e];
        this->a = a;
        this->e = e;
        this->lrh_l_ary = lrh_l_ary;
        this->lrh_l_map = lrh_l_map;
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
    virtual void solve() {
        throw "Should override solve()";
    }
    virtual void getSol() {
        throw "Should override getSol()";
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
              double ***lrh_l_ary, std::map<iiiTup, double> *lrh_l_map,
              std::string rutLogPath, std::string rutLpLogPath): RouterSCD(prob, tt, time_limit_sec, a, e, lrh_l_ary, lrh_l_map) {
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
    void solve();
    void getSol();
};

class RouterBnC: public RouterSCD {
public:
    const std::vector<std::string> cut_names {"hSE", "hCA", "hRS", "hIP"};
    const IloCplex::CutManagement cutManagerType = IloCplex::UseCutForce;
    const IloBool isLocalCutAdd = IloFalse;
    //
    rmm::BnC *rutAlgo;
    rgh::InsertionHeuristic* gh;
    //
    RouterBnC(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec,
              int a, int e,
              double ***lrh_l_ary, std::map<iiiTup, double> *lrh_l_map,
              std::string rutLogPath, std::string rutLpLogPath): RouterSCD(prob, tt, time_limit_sec, a, e, lrh_l_ary, lrh_l_map) {
        rutAlgo = new rmm::BnC(rutProb, tt, time_limit_sec, rutLogPath, rutLpLogPath, cut_names, true, cutManagerType, isLocalCutAdd);
        rutAlgo->cplex->setWarning(rutAlgo->env.getNullStream());
        this->rutCplex = rutAlgo->cplex;
        gh = new rgh::InsertionHeuristic(rutProb, tt, time_limit_sec, rutLogPath);
    }
    ~RouterBnC() {
        delete rutAlgo;
        delete gh;
    }
    //
    void update();
    void solve();
    void getSol();
};

class RouterBnCbc: public RouterBnC {
public:
    //
    RouterBnCbc(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec,
              int a, int e,
              double ***lrh_l_ary, std::map<iiiTup, double> *lrh_l_map,
              std::string rutLogPath, std::string rutLpLogPath): RouterBnC(prob, tt, time_limit_sec, a, e, lrh_l_ary, lrh_l_map, rutLogPath, rutLpLogPath) {
    }
    //
    void update();
};

class RouterGH: public RouterSCD {
public:
    rgh::InsertionHeuristic *rutAlgo;
    //
    RouterGH(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec,
            int a, int e,
            double ***lrh_l_ary, std::map<iiiTup, double> *lrh_l_map,
            std::string rutLogPath): RouterSCD(prob, tt, time_limit_sec, a, e, lrh_l_ary, lrh_l_map) {
        rutAlgo = new rgh::InsertionHeuristic(rutProb, tt, time_limit_sec, rutLogPath);
    }
    ~RouterGH() {
        delete rutAlgo;
    }
    //
    void update();
    void solve();
    void getSol();
};


#endif /* RouterSCD_h */
