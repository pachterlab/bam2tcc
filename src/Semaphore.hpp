#ifndef __SEMAPHORE_H__
#define __SEMAPHORE_H__

#include <condition_variable>
#include <mutex>

/*
 * Semaphore implementation adapted from
 * http://stackoverflow.com/questions/4792449/c0x-has-no-semaphores-how-to-synchronize-threads
 * by TA or prof of Wi 2018 CS2 at Caltech.
 */
class Semaphore {
private:
    std::mutex m;
    std::condition_variable cv;
    unsigned long count;

public:
    // Default to a mutex
    Semaphore() : count(1) {}

    Semaphore(unsigned long count) : count(count) {}

    ~Semaphore() {};

    void inc() {
        std::unique_lock<std::mutex> lock(m);
        count++;
        cv.notify_one();
    }

    void dec() {
        std::unique_lock<std::mutex> lock(m);
        while (!count)
            cv.wait(lock);
        count--;
    }
};

#endif
