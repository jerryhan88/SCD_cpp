//
//  ThreadPool.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 1/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ThreadPool.hpp"

void ThreadPool::resize(int newTCount) {
    int tmp = MAX_THREADS;
    if(newTCount > tmp || newTCount < 1){
        tmp = numThreads;
        numThreads = MAX_THREADS;
        Pool.resize(newTCount);
        for (int i = tmp; i != numThreads; ++i) {
            Pool.emplace_back(std::thread(&ThreadPool::threadManager, this));
            Pool.back().detach();
        }
    } else if (newTCount > numThreads) {
        uint8_t tmp = numThreads;
        numThreads = newTCount;
        Pool.resize(numThreads);
        for (int i = tmp; i != numThreads; ++i) {
            Pool.emplace_back(std::thread(&ThreadPool::threadManager, this));
            Pool.back().detach();
        }
    } else {
        numThreads = (uint8_t)newTCount;
        Pool.resize(newTCount);
    }
}

uint8_t ThreadPool::getThreadCount(){
    return numThreads;
}

void ThreadPool::threadManager() {
    while (true) {
        std::unique_lock<std::mutex> lock(JobMutex);
        thread_cond.wait(lock, [this] {return !JobQueue.empty(); });
        if (JobQueue.size() < 1)
            continue;
        std::shared_ptr<Job> job = JobQueue.front();
        JobQueue.pop();
        lock.release()->unlock();
        //
        job->execute();
        
    }
}
