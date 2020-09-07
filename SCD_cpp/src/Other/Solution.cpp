//
//  Solution.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 18/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Solution.hpp"

#define ROW_BUFFER_SIZE 2048


void Solution::writeSolCSV(std::string solPathCSV) {
    char row[ROW_BUFFER_SIZE];
    std::fstream fout_csv;
    fout_csv.open(solPathCSV, std::ios::out);
    fout_csv << "objV,Gap,elaCpuTime,elaWallTime,notes" << "\n";
    sprintf(row, "%f,%f,%f,%f,%s",
            objV, gap, cpuT, wallT, note.c_str());
    fout_csv << row << "\n";
    fout_csv.close();
}

void Solution::writeSolTXT(std::string solPathTXT) {
    char buf[2048];
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
        std::vector<std::pair<int, int>> associatedNodes;
        for (int k: prob->K) {
            if (y_ary[a][k] > 0.5) {
                assignedTasks.push_back(k);
                associatedNodes.push_back(std::pair<int, int> {prob->h_k[k], prob->n_k[k]});
            }
        }
        fout_txt << "A" << a << ": T [";
        for (int i = 0; i < assignedTasks.size(); i++) {
            if (i == 0) {
                fout_txt << assignedTasks[i];
            } else {
                fout_txt << ", " << assignedTasks[i];
            }
        }
        fout_txt << "]; N [";
        for (int i = 0; i < associatedNodes.size(); i++) {
            sprintf(buf, "(%d, %d)", associatedNodes[i].first, associatedNodes[i].second);
            if (i == 0) {
                fout_txt << buf;
            } else {
                fout_txt << ", " << buf;
            }
        }
        fout_txt << "]\n";
        //
        for (int e : prob->E_a[a]) {
            std::map<int, int> fromToPairs;
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    if (x_ary[a][e][i][j] > 0.5) {
                        fromToPairs[i] = j;
                    }
                }
            }
            int n0 = prob->o_ae[a][e];
            std::vector<int> route;
            while (n0 != prob->d_ae[a][e]) {
                route.push_back(n0);
                n0 = fromToPairs[n0];
            }
            route.push_back(n0);
            fout_txt << "\t R" << e << ": RR ";
            for (int i = 0; i < prob->R_ae[a][e].size(); i++) {
                if (i == 0) {
                    fout_txt << prob->R_ae[a][e][i];
                } else {
                    fout_txt << "-" << prob->R_ae[a][e][i];
                }
            }
            fout_txt << "; CT [";
            std::set<int> routeS(route.begin(), route.end());
            std::vector<int> deliveredTask;
            for (int k: assignedTasks) {
                if (routeS.find(prob->h_k[k]) != routeS.end() && routeS.find(prob->n_k[k]) != routeS.end()){
                    deliveredTask.push_back(k);
                }
            }
            for (int i = 0; i < deliveredTask.size(); i++) {
                if (i == 0) {
                    fout_txt << deliveredTask[i];
                } else {
                    fout_txt << ", " << deliveredTask[i];
                }
            }
            fout_txt << "]\n";
            fout_txt << "\t\t";
            for (int i = 0; i < route.size(); i++) {
                sprintf(buf, "%d(%.2f)", route[i], u_ary[a][e][route[i]]);
                if (i == 0) {
                    fout_txt << buf;
                } else {
                    fout_txt << "-" << buf;
                }
            }
            fout_txt << "\n";
        }
    }
    fout_txt.close();
}
