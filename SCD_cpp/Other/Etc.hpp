//
//  Etc.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 16/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef Etc_hpp
#define Etc_hpp

#include <dirent.h>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <chrono>

bool hasOption(std::vector<std::string> &arguments, std::string option);

std::string valueOf(std::vector<std::string> &arguments, std::string option);

std::vector<std::string> parseWithDelimiter(std::string str, std::string delimiter);

std::vector<std::string> read_directory(const std::string &d_path, const std::string &extension);


class TimeTracker {
public:
    std::clock_t c_start;
    std::chrono::high_resolution_clock::time_point w_start;
    //
    TimeTracker() {
        c_start = std::clock();
        w_start = std::chrono::high_resolution_clock::now();
    }
    ~TimeTracker() {}
    //
    std::string get_curTime();
    double get_elipsedTimeCPU();
    double get_elipsedTimeWall();
};

class FilePathOrganizer {
public:
    std::string logPath, solPathCSV, solPathTXT, lpPath, ilpPath;
    //
    FilePathOrganizer(const std::string &appr_dpath, const std::string &postfix);
    ~FilePathOrganizer() {}
};
#endif /* Etc_hpp */
