//
// Created by senbaikang on 18.05.21.
//

#ifndef DATAFILTER_THREAD_POOL_H
#define DATAFILTER_THREAD_POOL_H

#include <future>
#include <queue>
#include <thread>

/**
 * An abstract class of thread pool.
 */
class ThreadPool {

private:
  std::vector<std::thread> threads;
  std::queue<std::packaged_task<void()>> taskQueue;

  // Synchronization.
  std::mutex queueMutex;
  std::condition_variable job;
  std::condition_variable monitor;
  u_int32_t busy;
  u_int32_t processed;
  bool stop;

public:
  ThreadPool();
  explicit ThreadPool(const u_int16_t &);

  ~ThreadPool();

  template <class F, class... Args>
  decltype(auto) addTask(F &&, Args &&...);

  u_int32_t waitUntilFinished();

};

inline
ThreadPool::ThreadPool():
    ThreadPool(1)
{}

inline
ThreadPool::ThreadPool(const u_int16_t &t):
    busy(0),
    processed(0),
    stop(false)
{
  const u_int16_t providedThreadNum = t == 0 ? 1 : t;

  const u_int16_t _maxThreadNum = std::thread::hardware_concurrency();
  const u_int16_t maxThreadNum = _maxThreadNum == 0 ? 1 : _maxThreadNum;

  const u_int16_t &threadNum = std::min(providedThreadNum, maxThreadNum);

  for (u_int16_t i = 0; i < threadNum; i++)
  {
    threads.emplace_back(
        [this]()
        {
          for(;;)
          {
            std::packaged_task<void()> task;

            std::unique_lock<std::mutex> lock(queueMutex);

            job.wait(
                lock,
                [&]()
                {
                  return stop || !taskQueue.empty();
                }
                );

            if (stop && taskQueue.empty())
              break;

            ++busy;
            task = std::move(taskQueue.front());
            taskQueue.pop();
            lock.unlock();

            task();

            lock.lock();
            ++processed;
            --busy;
            monitor.notify_one();
          }
        }
        );
  }
}

inline
ThreadPool::~ThreadPool()
{
  {
    std::unique_lock<std::mutex> lock(queueMutex);
    stop = true;
  }

  job.notify_all();

  for (std::thread & i : threads)
    i.join();
}

template <class F, class... Args>
decltype(auto) ThreadPool::addTask(F &&f, Args &&...args)
{
  using return_type = decltype(f(args...));

  std::packaged_task<return_type()> task (
      std::bind(std::forward<F>(f), std::forward<Args>(args)...)
      );

  std::future<return_type> ret = task.get_future();

  {
    std::unique_lock<std::mutex> lock(queueMutex);

    if (stop)
      throw std::runtime_error("Error! Thread pool has been closed. No more tasks should be submitted.");

    taskQueue.emplace(std::move(task));
  }

  job.notify_one();

  return ret;
}

inline
u_int32_t ThreadPool::waitUntilFinished()
{
  std::unique_lock<std::mutex> lock(queueMutex);

  monitor.wait(
      lock,
      [this]() {
        return busy == 0 && taskQueue.empty();
      }
 );

  return processed;
}

#endif // DATAFILTER_THREAD_POOL_H
