//
//  main.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 16/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include <iostream>
#include <fstream>

#include <sys/types.h>


#include <vector>

#include "Other/Problem.hpp"
#include "Other/Solution.hpp"
#include "Other/Etc.hpp"
#include "Approach/Base.hpp"


int main(int argc, const char * argv[]) {
    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++) {
        arguments.push_back(argv[i]);
    }
    std::string basicFlags[] = {"-i"};
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
    //
    std::string prob_dpath(valueOf(arguments, "-i"));
    std::string appr_name("ILP");
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
        Problem *prob;
        prob = Problem::read_json(prob_fpath);
        std::string postfix = prob->problemName;
        FilePathOrganizer fpo(appr_dpath, postfix);
        if (!hasOption(arguments, "-l")) {
            fpo.logPath = "";
        }
        TimeTracker tt;
        std::cout << tt.get_curTime();
        std::cout << "\t" << appr_name << "; " << postfix << std::endl;
        //
        Base *solAppr;
        Solution *sol;
        if (appr_name_base == "ILP") {
            solAppr = new ILP(prob, fpo.logPath, &tt, fpo.lpPath);
            
            
            
            sol = solAppr->solve();
        } else {
            
        }
        
        sol->writeSolCSV(fpo.solPathCSV);
        sol->writeSolTXT(fpo.solPathTXT);
//        write_solution(prob, fpo, tt, mm);
    }
}
