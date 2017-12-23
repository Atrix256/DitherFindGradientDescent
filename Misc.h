#pragma once

#include <thread>
#include <vector>

//======================================================================================

inline float Lerp (float A, float B, float t)
{
    return A * (1.0f - t) + B * t;
}

//======================================================================================

template <typename T>
inline T Clamp (T A, T min, T max)
{
    if (A < min)
        return min;
    else if (A > max)
        return max;
    else return A;
}

//======================================================================================

template <typename THREADSTORAGE, typename LAMBDA>
void ForkJoin (std::vector<THREADSTORAGE>& threadStorage, const LAMBDA& lambda)
{
    // make an array of threads
    size_t numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(numThreads);

    // make storage for the threads
    threadStorage.resize(numThreads);

    // spin up the threads, calling the lambda
    for (size_t i = 0; i < numThreads; ++i)
    {
        threads[i] = std::thread(
            [i, &threadStorage, &lambda] ()
            {
                lambda(threadStorage[i]);
            }
        );
    }
 
    // wait for the threads to finish
    for (std::thread& t : threads)
        t.join();
}

//======================================================================================

template <typename LAMBDA>
void ForkJoin (const LAMBDA& lambda)
{
    // make an array of threads
    size_t numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(numThreads);

    // spin up the threads, calling the lambda
    for (size_t i = 0; i < numThreads; ++i)
    {
        threads[i] = std::thread(
            [&lambda] ()
            {
                lambda();
            }
        );
    }
 
    // wait for the threads to finish
    for (std::thread& t : threads)
        t.join();
}