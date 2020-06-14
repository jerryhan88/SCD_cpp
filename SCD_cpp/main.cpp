//
//  main.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 16/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>


#include <sys/types.h>
#include <limits>

#include <vector>

#include "include/Problem.hpp"
#include "include/Solution.hpp"
#include "include/SolApprBase.hpp"
#include "include/OtherSolvers.hpp"
#include "include/LRH.hpp"
#include "include/CGLC.hpp"
#include "include/ThreadPool.hpp"

#include "ck_util/util.hpp"         // from util
#include "ck_route/Other.hpp"     // from BnC_CPLEX


int main(int argc, const char * argv[]) {
    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++) {
        arguments.push_back(argv[i]);
    }
    std::string basicFlags[] = {"-i", "-a"};
    int numBasicFlags = sizeof(basicFlags) / sizeof(std::string);
    for (int i = 0; i < numBasicFlags; i++) {
        std::string option = basicFlags[i];
        if (!hasOption(arguments, option)) {
            std::string msg("Please provide the basic arguments;");
            msg += " " + option;
            std::cout << msg << std::endl;
            return 1;
        }
    }
    int numThreads = 1;
    unsigned long time_limit_sec = ULONG_MAX;
    std::string lp_algo = "DSM";
    bool enforcementMode = false;
    if (hasOption(arguments, "-nth")) {
        numThreads = std::stoi(valueOf(arguments, "-nth"));
    }
    if (hasOption(arguments, "-t")) {
        time_limit_sec = std::stoi(valueOf(arguments, "-t"));
    }
    if (hasOption(arguments, "-la")) {
        lp_algo = std::stoi(valueOf(arguments, "-la"));
    }
    if (hasOption(arguments, "-ef")) {
        enforcementMode = true;
    }
    //
    std::string prob_dpath(valueOf(arguments, "-i"));
    std::string appr_name(valueOf(arguments, "-a"));
    std::string appr_dpath = prob_dpath.substr(0, prob_dpath.find_last_of("/"));
    appr_dpath += "/" + appr_name;
    system(("mkdir -p " + appr_dpath).c_str());
    //
    std::string appr_name_base = appr_name.substr(0, appr_name.find_first_of("-"));
    std::vector<std::string> probFileNames = read_directory(prob_dpath, ".json");
    for (std::string fn: probFileNames) {
        std::string prob_fpath(prob_dpath + "/" + fn);
        std::cout << prob_fpath << std::endl;
        //
        std::cout << "\tprob_read_start; " << TimeTracker::get_curTime();
        Problem *prob;
        prob = Problem::read_json(prob_fpath);
        prob->gen_aeProbs();
        std::string postfix = prob->problemName;
        FilePathOrganizer fpo(appr_dpath, postfix);
        std::cout << "\tprob_read_end; " << TimeTracker::get_curTime();
        if (!enforcementMode) {
            std::ifstream is;
            std::string handledFiles[] = {fpo.solPathCSV, fpo.lpPath};
            int numAssoFiles = sizeof(handledFiles) / sizeof(std::string);
            bool isHandled = false;
            for (int i = 0; i < numAssoFiles ; i++) {
                is.open(handledFiles[i]);
                if (!is.fail()) {
                    is.close();
                    isHandled = true;
                    break;
                }
            }
            if (isHandled) {
                continue;
            }
        }
        if (!hasOption(arguments, "-l")) {
            fpo.logPath = "";
        }
        std::cout << "\t" << appr_name << "; " << postfix << "; " << TimeTracker::get_curTime();
        //
        TimeTracker tt;
        SolApprBase *solAppr = nullptr;
        Solution *sol;
        if (appr_name_base == "LP") {
            solAppr = new LP(prob, &tt, time_limit_sec, numThreads, fpo.logPath, fpo.lpPath, lp_algo);
            sol = solAppr->solve();
            sol->writeSolCSV(fpo.solPathCSV);
            delete solAppr; delete sol;
            continue;
        } else if (appr_name_base == "ILP") {
            solAppr = new ILP(prob, &tt, time_limit_sec, numThreads, fpo.logPath, fpo.lpPath, lp_algo);
            sol = solAppr->solve();
        } else if (appr_name_base == "LRH" || appr_name_base == "CGLC") {
            std::vector<std::string> tokens = parseWithDelimiter(appr_name, "-");
            assert (2 == tokens.size());
            std::string _router = tokens[1];
            char buf[4096];
            std::ifstream is;
            std::string setting_fpath(appr_dpath + "/LRH_setting.txt");
            is.open(setting_fpath);
            if (!is.fail()) {
                sprintf(buf, "rm %s", setting_fpath.c_str());
                system(buf);
            }
            double dual_gap_limit;
            unsigned int num_iter_limit, no_improvement_limit;
            if (hasOption(arguments, "-du")) {
                dual_gap_limit = std::stod(valueOf(arguments, "-du"));
            } else {
                dual_gap_limit = 0.05;
            }
            sprintf(buf, "echo DUAL_GAP_LIMIT: %f >> %s", dual_gap_limit, setting_fpath.c_str());
            system(buf);
            //
            if (hasOption(arguments, "-ni")) {
                num_iter_limit = std::stod(valueOf(arguments, "-ni"));
            } else {
                num_iter_limit  = 20;
            }
            sprintf(buf, "echo NUM_ITER_LIMIT: %d >> %s", num_iter_limit, setting_fpath.c_str());
            system(buf);
            //
            if (hasOption(arguments, "-im")) {
                no_improvement_limit = std::stoi(valueOf(arguments, "-im"));
            } else {
                no_improvement_limit = 3;
            }
            sprintf(buf, "echo NO_IMPROVEMENT_LIMIT: %d >> %s", no_improvement_limit, setting_fpath.c_str());
            system(buf);
            sprintf(buf, "echo TIME_LIMIT_SEC: %lu >> %s", time_limit_sec, setting_fpath.c_str());
            system(buf);
            //
            if (appr_name_base == "LRH") {
                solAppr = new LRH(prob, &tt,
                time_limit_sec, numThreads,
                fpo.logPath, fpo.lpPath, lp_algo,
                _router,
                dual_gap_limit, num_iter_limit, no_improvement_limit);
            } else {
                assert(appr_name_base == "CGLC");
                solAppr = new CGLC(prob, &tt,
                                time_limit_sec, numThreads,
                                fpo.logPath, fpo.lpPath, lp_algo,
                                _router,
                                dual_gap_limit, num_iter_limit, no_improvement_limit);
            }
            sol = solAppr->solve();
        } else if (appr_name_base == "PureGH") {
            solAppr = new PureGH(prob, &tt, time_limit_sec, numThreads, fpo.logPath);
            
            sol = solAppr->solve();
        } else {
            sol = nullptr;
        }
        sol->writeSolCSV(fpo.solPathCSV);
        sol->writeSolTXT(fpo.solPathTXT);
        delete solAppr;
        delete sol;
    }
}
