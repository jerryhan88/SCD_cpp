//
//  ThreadPool.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 1/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef ThreadPool_hpp
#define ThreadPool_hpp

#include <stdio.h>

#endif /* ThreadPool_hpp */



// https://codereview.stackexchange.com/questions/221617/thread-pool-c-implementation

#include <iostream>

#pragma once

#include<thread>
#include<vector>
#include<queue>
#include<mutex>
#include<condition_variable>
#include<functional>
#include<future>

#define MAX_THREADS std::thread::hardware_concurrency() - 1;


class Job {
public:
    virtual ~Job() {}
    virtual void execute() = 0;
};

template <typename RetType>
class AnyJob : public Job {
private:
    std::packaged_task<RetType()> func;
public:
    AnyJob(std::packaged_task<RetType()> func) : func(std::move(func)) {}
    void execute() {
        func();
    }
};

class ThreadPool{
public:
    static ThreadPool& getInstance(int numThreads) {
        static ThreadPool instance(numThreads);
        return instance;
    }
    template <typename Func, typename... Args >
    auto push(Func&& f, Args&&... args){
        typedef decltype(f(args...)) retType;
        std::packaged_task<retType()> task(std::move(std::bind(f, args...)));
        std::unique_lock<std::mutex> lock(JobMutex);
        std::future<retType> future = task.get_future();
        JobQueue.emplace( std::make_shared<AnyJob<retType> > (std::move(task)) );
        thread_cond.notify_one();
        return future;
    }
    //
    void resize(int newTCount);
    uint8_t getThreadCount();
    //
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool &operator=(const ThreadPool&) = delete;
private:
    uint8_t numThreads; // number of threads in the pool
    std::vector<std::thread> Pool; //the actual thread pool
    std::queue<std::shared_ptr<Job>> JobQueue;
    std::condition_variable thread_cond;// used to notify threads about available jobs
    std::mutex JobMutex; // used to push/pop jobs to/from the queue
    //
    void threadManager();
    //
    ThreadPool(uint8_t numThreads) : numThreads(numThreads) {
        int tmp = MAX_THREADS;
        if(numThreads > tmp){
            numThreads = tmp;
        }
        Pool.reserve(numThreads);
        for(int i = 0; i != numThreads; ++i){
            Pool.emplace_back(std::thread(&ThreadPool::threadManager, this));
            Pool.back().detach();
        }
    }
};

