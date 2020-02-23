//
//  Solution.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 18/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "Solution.hpp"

#define ROW_BUFFER_SIZE 2048


Solution::Solution(Problem *prob) {
    this-> prob = prob;
}

Solution::~Solution() {
    size_t numNodes = prob->locID.size();
    for (int a: prob->A) {
        delete [] y_ak[a];
        for (int e: prob->E_a[a]) {
            delete [] z_aek[a][e];
            delete [] u_aei[a][e];
            for (int i = 0; i < numNodes; i++) {
                delete [] x_aeij[a][e][i];
            }
            delete [] x_aeij[a][e];
        }
        delete [] z_aek[a];
        delete [] u_aei[a];
        delete [] x_aeij[a];
    }
    delete [] y_ak;
    delete [] z_aek;
    delete [] x_aeij;
    delete [] u_aei;
}

void Solution::alloMem4dvs() {
    size_t numAgents, numTasks, numNodes;
    numAgents = prob->A.size();
    numTasks = prob->K.size();
    numNodes = prob->locID.size();
    //
    y_ak = new double*[numAgents];
    z_aek = new double**[numAgents];
    x_aeij = new double***[numAgents];
    u_aei = new double**[numAgents];
    for (int a: prob->A) {
        y_ak[a] = new double[numTasks];
        size_t numRR = prob->E_a[a].size();
        z_aek[a] = new double*[numRR];
        x_aeij[a] = new double**[numRR];
        u_aei[a] = new double*[numRR];
        for (int e: prob->E_a[a]) {
            z_aek[a][e] = new double[numTasks];
            x_aeij[a][e] = new double*[numNodes];
            u_aei[a][e] = new double[numNodes];
            for (int i = 0; i < numNodes; i++) {
                x_aeij[a][e][i] = new double[numNodes];
            }
        }
    }
}

void Solution::writeSolCSV(std::string solPathCSV) {
    char row[ROW_BUFFER_SIZE];
    std::fstream fout_csv;
    fout_csv.open(solPathCSV, std::ios::out);
    fout_csv << "objV,Gap,eliCpuTime,eliWallTime,notes" << "\n";
    sprintf(row, "%f,%f,%f,%f,%s",
            objV, gap, cpuT, wallT, note.c_str());
    fout_csv << row << "\n";
    fout_csv.close();
}

void Solution::writeSolTXT(std::string solPathTXT) {
    std::fstream fout_txt;
    fout_txt.open(solPathTXT, std::ios::out);
    //
    fout_txt << "Summary\n";
    fout_txt << "\t ObjV: " << objV << "\n";
    fout_txt << "\t Gap: " << gap << "\n";
    fout_txt << "\t Cpu Time: " << cpuT << "\n";
    fout_txt << "\t Wall Time: " << wallT << "\n";
    fout_txt << "\n";
    //
    fout_txt << "Details\n";
    for (int a: prob->A) {
        std::vector<int> assignedTasks;
        std::vector<std::string> _assignedTasks;
        std::set<int> meaninglessNodes;
        for (int i: prob->N) {
            meaninglessNodes.insert(i);
        }
        for (int k: prob->K) {
            if (y_ak[a][k] > 0.5) {
                assignedTasks.push_back(k);
                _assignedTasks.push_back(std::to_string(k));
                if (meaninglessNodes.find(prob->h_k[k]) != meaninglessNodes.end()) {
                    meaninglessNodes.erase(prob->h_k[k]);
                }
                if (meaninglessNodes.find(prob->n_k[k]) != meaninglessNodes.end()) {
                    meaninglessNodes.erase(prob->n_k[k]);
                }
            }
        }
        fout_txt << "A" << a << ": [";
        for (int i = 0; i < _assignedTasks.size(); i++) {
            if (i != 0) {
                fout_txt << ",";
            }
            fout_txt << _assignedTasks[i];
        }
        fout_txt << "]\n";
        char buf[1024];
        for (int e : prob->E_a[a]) {
            sprintf(buf, "s0_%d_%d", a, e);
            int ae_o = prob->locID[buf];
            sprintf(buf, "s%d_%d_%d", (int) prob->S_ae[a][e].size(), a, e);
            int ae_d = prob->locID[buf];
            std::map<int, int> _route;
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    if (x_aeij[a][e][i][j] > 0.5) {
                        _route[i] = j;
                    }
                }
            }
            int i = ae_o;
            std::ostringstream route;
            std::vector<std::string> _accomplishedTasks;
            while (i != ae_d) {
                if (meaninglessNodes.find(i) == meaninglessNodes.end()) {
                    sprintf(buf, "%s(%.2f)-", prob->idLoc[i].c_str(), u_aei[a][e][i]);
                    route << buf;
                    if (prob->idLoc[i].rfind("n", 0) == 0) {
                        _accomplishedTasks.push_back(prob->idLoc[i].substr(1));
                    }
                }
                i = _route[i];
            }
            sprintf(buf, "%s(%.2f)", prob->idLoc[i].c_str(), u_aei[a][e][i]);
            route << buf;
            fout_txt << "\t R" << e << "[";
            for (int i = 0; i < _accomplishedTasks.size(); i++) {
                if (i != 0) {
                    fout_txt << ",";
                }
                fout_txt << _accomplishedTasks[i];
            }
            fout_txt << "]: " << route.str() << "\n";
        }
    }
    fout_txt.close();
}
