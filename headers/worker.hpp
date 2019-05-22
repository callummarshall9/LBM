#include <thread>
#include <mutex>
#include <chrono>

class worker
{
public:
  worker();

  // Thread function.
  void do_work(void (*notify)(float progress, std::chrono::milliseconds milliseconds, int runs_left));

  void get_data(double* fraction_done) const;
  void stop_work();
  bool has_stopped() const;

private:
  mutable std::mutex thread_mutex;
  bool finished;
  float fraction_done;
};
