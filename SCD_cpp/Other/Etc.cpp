//
//  Etc.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 16/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "Etc.hpp"


bool hasOption(std::vector<std::string> &arguments, std::string option) {
    bool hasValue = false;
    for (std::string str: arguments) {
        if (str == option) {
            hasValue = true;
        }
    }
    return hasValue;
}

std::string valueOf(std::vector<std::string> &arguments, std::string option) {
    std::string value = "";
    for (int i = 0; i < arguments.size(); i++) {
        if (arguments[i] == option) {
            value = arguments[i + 1];
        }
    }
    return value;
}

std::vector<std::string> parseWithDelimiter(std::string str, std::string delimiter) {
    std::vector<std::string> tokens;
    size_t pos = 0;
    std::string token;
    while ((pos = str.find(delimiter)) != std::string::npos) {
        token = str.substr(0, pos);
        tokens.push_back(token);
        str.erase(0, pos + delimiter.length());
    }
    tokens.push_back(str);
    return tokens;
    
}

std::vector<std::string> read_directory(const std::string &d_path, const std::string &extension) {
    std::vector<std::string> fileNames;
    DIR* dirp = opendir(d_path.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        std::string fn(dp->d_name);
        std::size_t found = fn.find(extension);
        if (found!=std::string::npos)
            fileNames.push_back(fn);
    }
    closedir(dirp);
    std::sort(fileNames.begin(), fileNames.end());
    return fileNames;
}

std::string TimeTracker::get_curTime() {
    time_t now = time(0);
    char* dt = ctime(&now);
    return std::string(dt);
}

double TimeTracker::get_elipsedTimeCPU() {
    std::clock_t c_end = std::clock();
    return (c_end-c_start) / (double) CLOCKS_PER_SEC;
}

double TimeTracker::get_elipsedTimeWall() {
    std::chrono::high_resolution_clock::time_point w_end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(w_end-w_start).count() / 1000.0;
}

FilePathOrganizer::FilePathOrganizer(const std::string &appr_dpath, const std::string &postfix) {
    std::string appr_name = appr_dpath.substr(appr_dpath.find_last_of("/") + 1, appr_dpath.size());
    std::string prefix(appr_dpath + "/" + appr_name + "_" + postfix);
    //
    this->logPath = prefix + ".log";
    this->solPathCSV = prefix + ".csv";
    this->solPathTXT = prefix + ".txt";
    this->lpPath = prefix + ".lp";
    this->ilpPath = prefix + ".ilp";
}

