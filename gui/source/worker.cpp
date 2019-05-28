#include "worker.hpp"
#include <sstream>
#include <chrono>
#include "headers/LBM.hpp"

worker::worker() :
  thread_mutex(),
  finished(false),
  fraction_done(0.0)
{
}

// Accesses to these data are synchronized by a mutex.
// Some microseconds can be saved by getting all data at once, instead of having
// separate get_fraction_done() and get_message() methods.
void worker::get_data(double* fraction_done_update) const
{
  std::lock_guard<std::mutex> lock(thread_mutex);

  if (fraction_done_update) {
    *fraction_done_update = fraction_done;
  }
}

void worker::stop_work()
{
  std::lock_guard<std::mutex> lock(thread_mutex);
  finished = true;
}

bool worker::has_stopped() const
{
  std::lock_guard<std::mutex> lock(thread_mutex);
  return finished;
}

void worker::do_work(void (*notify)(float fraction_done, std::chrono::milliseconds milliseconds, int runs_left))
{
  {
    std::lock_guard<std::mutex> lock(thread_mutex);
    finished = false;
    fraction_done = 0.0;
     // The mutex is unlocked here by lock's destructor.
  }
  LBM solver(64);
 	int scale = 1;
 	int runs = 200 * scale * scale * scale;
  auto last_call = std::chrono::high_resolution_clock::now();
 	for(int i = 0; i < runs; i = i + 1) {
    if(finished) { continue; }
    std::lock_guard<std::mutex> lock(thread_mutex);
 		fraction_done = (float)(i+1) / (float)(runs);
 		solver.output_lbm_data("output/" + std::to_string(i) + ".csv");
 		solver.perform_timestep();
    auto after_call = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(after_call - last_call);
    last_call = after_call;
    (*notify)(fraction_done, duration, runs - (i+1));
 	}
  {
    std::lock_guard<std::mutex> lock(thread_mutex);
    finished = true;
    std::chrono::milliseconds time(0);
    (*notify)(1.0, time, 0);
  }
}
