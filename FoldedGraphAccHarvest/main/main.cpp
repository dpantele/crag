#include "acc.h"

#include <bitset>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <future>
#include <list>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <queue>

using namespace crag;

std::vector<std::string> Split(const std::string& str, char split = ':') {
  std::vector<std::string> result;
  decltype(str.find(split)) last_split = 0;
  while (last_split < str.size()) {
    auto new_split = str.find(split, last_split);
    result.push_back(str.substr(last_split, new_split - last_split));
    if (new_split == std::string::npos) {
      break;
    }
    last_split = new_split + 1;
  }
  return result;
}

struct Stopwatch {
  using clock = std::chrono::high_resolution_clock;
  clock::duration duration_{0};
  clock::duration last_duration_{0};
  uint64_t iterations_count_ = 0u;

  typedef std::chrono::microseconds ResultUnits;

  struct Iteration {
    Stopwatch* stopwatch_;
    clock::duration duration_{0};
    clock::time_point last_click_;
    uint64_t clicks_count_{0};

    Iteration(Stopwatch* stopwatch)
      : stopwatch_(stopwatch)
    { }

    Iteration(const Iteration& other) = delete;

    Iteration(Iteration&& other)
      : stopwatch_(other.stopwatch_)
      , duration_(other.duration_)
      , last_click_(other.last_click_)
      , clicks_count_(other.clicks_count_)
    { 
      other.clicks_count_ = 0;
      other.duration_ = other.duration_.zero();
    }

    inline bool IsRunning() const {
      return clicks_count_ % 2 == 1;
    }

    void Click() {
      if (IsRunning()) {
        duration_ += (clock::now() - last_click_);
      } else {
        last_click_ = clock::now();
      }

      ++clicks_count_;
    }

    void Reset() {
      if (IsRunning()) {
        Click();
      }
      if (clicks_count_) {
        stopwatch_->Report(this);
        clicks_count_ = 0;
        duration_ = clock::duration(0);
      }
    }

    ~Iteration() {
      Reset();
    }
  };

  void Report(Iteration* iteration) {
    ++iterations_count_;
    duration_ += iteration->duration_;
    last_duration_ = iteration->duration_;
  }

  Iteration NewIter() {
    return Iteration(this);
  }

  long long last() {
    return std::chrono::duration_cast<ResultUnits>(last_duration_).count();
  }

  long long average() {
    return iterations_count_ ? std::chrono::duration_cast<ResultUnits>(duration_ / iterations_count_).count() : 0;
  }

  uint64_t iterations() const {
    return iterations_count_;
  }

  long long total_time() const {
    return duration_.count();
  }
};

template<typename U, typename V>
std::pair<V, U> Swapped(std::pair<U, V> p) {
  return std::pair<V, U>(std::move(p.second), std::move(p.first));
}

//structures and algoritms for main below:

struct Stopwatches {
  typedef Stopwatch Field;
  Field folding_;
  Field harvest_;
  Field normalize_;
  Field reweight_;

  Field& folding() {
    return folding_;
  }
  Field& harvest() {
    return harvest_;
  }
  Field& normalize() {
    return normalize_;
  }
  Field& reweight() {
    return reweight_;
  }

  Stopwatches() { }

  //non-movable, non-copyable
  Stopwatches(const Stopwatches&) = delete;
  Stopwatches(Stopwatches&&) = delete;
  Stopwatches& operator=(const Stopwatches&) = delete;
  Stopwatches& operator=(Stopwatches&&) = delete;

  long long total() const {
    long long total = 0;
    total += folding_.total_time();
    total += harvest_.total_time();
    total += normalize_.total_time();
    total += reweight_.total_time();
    return total;
  }
};

struct StopwatchesIteration {
  typedef Stopwatch::Iteration Field;
  Field folding_;
  Field harvest_;
  Field normalize_;
  Field reweight_;

  Field& folding() {
    return folding_;
  }
  Field& harvest() {
    return harvest_;
  }
  Field& normalize() {
    return normalize_;
  }
  Field& reweight() {
    return reweight_;
  }

  StopwatchesIteration(Stopwatches* s)
    : folding_(s->folding().NewIter())
    , harvest_(s->harvest().NewIter())
    , normalize_(s->normalize().NewIter())
    , reweight_(s->reweight().NewIter())
  { }

  StopwatchesIteration(const StopwatchesIteration&) = delete;
  StopwatchesIteration(StopwatchesIteration&&) = delete;
  StopwatchesIteration& operator=(const StopwatchesIteration&) = delete;
  StopwatchesIteration& operator=(StopwatchesIteration&&) = delete;

  void Reset() {
    folding_.Reset();
    harvest_.Reset();
    normalize_.Reset();
    reweight_.Reset();
  }
};

struct Parameters {
  size_t max_harvest_length = Word::kMaxLength;
  size_t workers_count_ = 1;
  std::vector<uint16_t> complete_count;
  std::pair<std::string, std::string> initial_strings;
  std::pair<std::string, std::string> required_strings;

  Parameters()
    : complete_count(Word::kMaxLength + 1, 2)
    , initial_strings("xyxYXY", "xxxYYYY")
    , required_strings("x", "y")
  { }
};

class PairToProcess {
 public:
  PairToProcess(Word u, Word v, const Parameters& p, Stopwatches* timer) 
    : task_id_(GetNextId())
    , u_(std::move(u))
    , v_(std::move(v))
    , p_(p)
    , time_(timer)
  { }

  PairToProcess(const PairToProcess&) = delete;
  PairToProcess(PairToProcess&&) = delete;
  PairToProcess& operator=(const PairToProcess&) = delete;
  PairToProcess& operator=(PairToProcess&&) = delete;
  
  const Word& u() const {
    return u_;
  }

  const Word& v() const {
    return v_;
  }

  int task_id() const {
    return task_id_;
  }

  struct Stats {
    int worker_id_ = 0;
    uint64_t graph_size_ = 0;
    uint64_t edges_before_reweight_ = 0;
    uint64_t edges_after_reweight_ = 0;
    uint64_t count_after_harvest_ = 0;
    uint64_t count_after_normalize_ = 0;
  };

  Stats s_;
  StopwatchesIteration time_;
  std::vector<std::pair<Word, Word>> generated_pairs_;

 
  void InWorkBy(int worker_id) {
    s_.worker_id_ = worker_id;
  }

  void Complete() {
    std::unique_lock<std::mutex> lk(m_computed_);
    computed_ = true;
    lk.unlock();
    cv_computed_.notify_one();
  }

  void Wait() {
    std::unique_lock<std::mutex> lk(m_computed_);
    while (!computed_) {
      cv_computed_.wait(lk);
    }
  }

  void Process() {
    time_.folding().Click();
    FoldedGraph2 g;
    g.PushCycle(u(), g.root(), 1);

    for (auto i = 0u; i < p_.complete_count[v().size()]; ++i) {
      g.CompleteWith(v());
    }
    time_.folding().Click();

    s_.graph_size_ = g.size();
    s_.edges_before_reweight_ = g.CountNontrivialEdges();

    time_.reweight().Click();
    g.Reweight();
    time_.reweight().Click();

    s_.edges_after_reweight_ = g.CountNontrivialEdges();

    time_.harvest().Click();
    auto eq_u = g.Harvest(p_.max_harvest_length, g.root());
    time_.harvest().Click();

    s_.count_after_harvest_ = eq_u.size();

    for (auto u_p = eq_u.begin(); u_p != eq_u.end(); ++u_p) {
      time_.normalize().Click();
      generated_pairs_.push_back(GetCanonicalPair(v(), *u_p)); 
      time_.normalize().Click();
    }
    std::sort(generated_pairs_.begin(), generated_pairs_.end());
    auto new_end = std::unique(generated_pairs_.begin(), generated_pairs_.end());
    generated_pairs_.erase(new_end, generated_pairs_.end());  
  }
 private:
  static int GetNextId() {
    static int next_task_id = 0;
    return ++next_task_id;
  }
  int task_id_; //!< Number in the queue
  Word u_;      //!< Generator
  Word v_;      //!< Relator

  const Parameters& p_; //!< General parameters of the algorithm

  bool computed_ = false;
  std::condition_variable cv_computed_;
  mutable std::mutex m_computed_;
};

template<typename T>
class shared_queue
{
  std::queue<T> queue_;
  mutable std::mutex m_;
  std::condition_variable data_cond_;

  shared_queue& operator=(const shared_queue&) = delete;
  shared_queue(const shared_queue& other) = delete;

public:
  shared_queue(){}

  void push(T item){
    std::lock_guard<std::mutex> lock(m_);
    queue_.push(std::move(item));
    data_cond_.notify_one();
  }

  /// \return immediately, with true if successful retrieval
  bool try_and_pop(T& popped_item){
    std::lock_guard<std::mutex> lock(m_);
    if(queue_.empty()){
      return false;
    }
    popped_item = std::move(queue_.front());
    queue_.pop();
    return true;
  }

  /// Try to retrieve, if no items, wait till an item is available and try again
  bool wait_and_pop(T& popped_item, std::chrono::milliseconds timeout){
    std::unique_lock<std::mutex> lock(m_); // note: unique_lock is needed for std::condition_variable::wait
    if (data_cond_.wait_for(lock, timeout, [this](){ return !queue_.empty(); })) {
      popped_item = std::move(queue_.front());
      queue_.pop();
      return true;
    }
    return false;
  }

  bool empty() const{
    std::lock_guard<std::mutex> lock(m_);
    return queue_.empty();
  }

  unsigned size() const{
    std::lock_guard<std::mutex> lock(m_);
    return queue_.size();
  }
};
template<typename Task>
class Worker {
public:
  shared_queue<Task> *const task_;
  bool exit_ = false;

  std::thread my_thread_;
  int worker_id_;
 private:
  Worker(int worker_id, shared_queue<Task> *const task)
    : task_(task)
    , worker_id_(worker_id)
  { }

public:
  static std::unique_ptr<Worker> Create(shared_queue<Task> *const task) {
    static int worker_id = 0;
    std::unique_ptr<Worker> worker(new Worker(++worker_id, task));
    worker->my_thread_ = std::thread(&Worker::Run, worker.get());
    return worker;
  }
  
  Worker(const Worker&) = delete;

  void Run() {
    while (!exit_) {
      Task t = nullptr;
      if (!task_->wait_and_pop(t, std::chrono::milliseconds{100})) {
        continue;
      }
      assert(t);

      t->InWorkBy(worker_id_);
      t->Process();
      t->Complete();
    }
  }

  ~Worker() {
    exit_ = true;
    my_thread_.join();
  }
};

int main(int argc, const char *argv[]) {
  Parameters p;
  std::ostream* out = &std::cout;

  std::ofstream stats_out;
  std::ofstream unproc_words;
  std::ofstream proc_words;
  std::ofstream estats_out;

  std::ifstream words_in("h15_no_auto_ak3_unproc_words.txt");

  std::pair<std::string, std::string> next_pair;

  words_in.ignore(64, '\n');

  std::set<std::pair<Word, Word>> though_for_words;

  while (words_in) {
    std::string appeared_in;
    words_in >> appeared_in;
    words_in >> next_pair.first;
    next_pair.first.pop_back();
    words_in >> next_pair.second;
    if (words_in) {
      though_for_words.emplace(Word(next_pair.first), Word(next_pair.second));
    }
  }


  for (int argi = 1; argi < argc; ++argi) {
    auto arg = Split(argv[argi]);
    if (arg.empty()) {
      continue;
    }

    if (arg.front() == "maxhl") {
      p.max_harvest_length = std::stoi(arg[1]);
    } else if (arg.front() == "j") {
      p.workers_count_ = std::stoi(arg[1]);
    } else if (arg.front() == "comp") {
      auto up_to = std::stoi(arg[1]);
      auto count = std::stoi(arg[2]);
      if (p.complete_count.size() < up_to) {
        p.complete_count.resize(up_to, count);
      }
      for (size_t i = 0; i < up_to; ++i) {
        p.complete_count[i] = count; 
      }
    } else if (arg.front() == "init") {
      p.initial_strings.first = arg[1];
      p.initial_strings.second = arg[2];
    } else if (arg.front() == "out") {
      stats_out.open(arg[1] + ".txt");
      out = &stats_out;

      proc_words.open(arg[1] + "_proc_words.txt");
      proc_words << "iter, u, v\n";
      unproc_words.open(arg[1] + "_unproc_words.txt");
      unproc_words << "u, v\n";
      estats_out.open(arg[1] + "_ext.txt");
    }
  }

  *out << std::left ;
  *out << std::setw(7) << "iter";
  *out << std::setw(9) << ", unproc";
  *out << std::setw(9) << ", total";
  *out << std::setw(4) << ", u";
  *out << std::setw(4) << ", v";
  *out <<                 ", added pairs";
  *out <<                 ", folding time";
  *out <<                 ", harvest time";
  *out << std::endl;

  if (estats_out.is_open()) {
    estats_out << "Configuration: " << std::endl;
    estats_out << "Max length for harvest: " << p.max_harvest_length << std::endl;
    estats_out << "Pairs are also normalized using automorphisms: yes";
    estats_out << "Initial words:          " << p.initial_strings.first << " | " << p.initial_strings.second << std::endl;
    estats_out << "\nHow many times graph is completed with v: " << std::endl;
    estats_out << "Length of v: ";
    for (auto i = 0u; i < p.complete_count.size(); ++i) {
      estats_out << std::setw(2) << i << ' ';
    }
    estats_out << "\nTimes:       ";
    for (auto i = 0u; i < p.complete_count.size(); ++i) {
      estats_out << std::setw(2) << p.complete_count[i] << ' ';
    }
    estats_out << "\n\n";
    estats_out << "iteration, unprocessed count, total count, ";
    estats_out << "graph size, ";
    estats_out << "non-zero edges before reweight, ";
    estats_out << "non-zero edges after reweight, ";
    estats_out << "harvest result size, ";
    estats_out << "after normalize size, ";
    estats_out << "pairs added count, ";
    estats_out << "u size, v size, ";
    estats_out << "worker id, ";
    estats_out << "folding time, average folding time, ";
    estats_out << "reweight time, average reweight time, ";
    estats_out << "harvest time, average harvest time, ";
    estats_out << "normalize time, average normalize time";
    estats_out << "\n";
  }

  auto initial = GetCanonicalPair(p.initial_strings.first.c_str(), p.initial_strings.second.c_str());
  auto required = GetCanonicalPair(p.required_strings.first.c_str(), p.required_strings.second.c_str());

  std::set<std::pair<Word, Word>> unprocessed_pairs = {};
  std::set<std::pair<Word, Word>> all_pairs = {};

  int counter = 0;

  shared_queue<PairToProcess*> tasks;
  Stopwatches total_time;

  std::vector<std::unique_ptr<Worker<PairToProcess*>>> workers;
  for (auto i = 0; i < p.workers_count_; ++i) {
    workers.push_back(Worker<PairToProcess*>::Create(&tasks));
  }

  Stopwatch overall_time;
  auto overall_time_stamp = overall_time.NewIter();
  overall_time_stamp.Click();

  std::deque<PairToProcess> future_results;
  auto NewTask = [&future_results, &unproc_words, &tasks, &p, &total_time](const std::pair<Word, Word>& pair, int placed_by) -> int {
    future_results.emplace_back(pair.first, pair.second, p, &total_time);
    tasks.push(&future_results.back());
    if (unproc_words.is_open()) {
      unproc_words << placed_by << ", ";
      PrintWord(pair.first, &unproc_words);
      unproc_words << ", ";
      PrintWord(pair.second, &unproc_words);
      unproc_words << "\n";
    }
    return future_results.back().task_id();
  };

  auto AddPair = [&NewTask, &all_pairs, &unproc_words, &counter, &though_for_words](const std::pair<Word, Word>& pair, int added_by) -> int {
    auto exists = all_pairs.insert(pair);
    if (exists.second) {
      if (though_for_words.count(pair)) {
        std::cout << "FOUND!!!!" << std::endl;
        PrintWord(pair.first, &std::cout);
        std::cout << std::endl;
        PrintWord(pair.second, &std::cout);
        std::cout << std::endl;

        exit(0);
      }

      return NewTask(pair, added_by);
    }
    return false;
  };

  auto initial_task_id = AddPair(initial, 0);
  AddPair(Swapped(initial), initial_task_id);

  while (!all_pairs.count(required) && !future_results.empty()) {
    future_results.front().Wait();

    auto& next_result = future_results.front();

    auto pairs_count_before = all_pairs.size();

    for (auto&& new_pair : next_result.generated_pairs_) {
      auto new_id = AddPair(new_pair, next_result.task_id());
      if (new_id) {
        AddPair(Swapped(new_pair), new_id);
      }
    }

    auto pairs_added = all_pairs.size() - pairs_count_before;
    next_result.time_.Reset();

    *out << std::left << std::setw(7) << next_result.task_id() << ", ";
    *out << std::right << std::setw(7) << pairs_count_before - next_result.task_id() << ", ";
    *out << std::right << std::setw(7) << pairs_count_before << ", ";
    *out << std::setw(2) << next_result.u().size() << ", ";
    *out << std::setw(2) << next_result.v().size() << ", ";
    *out << pairs_added << ", ";
    *out << total_time.folding().last() << ", ";
    *out << total_time.harvest().last() << std::endl;

    if (estats_out.is_open()) {
      estats_out << next_result.task_id() << ", ";
      estats_out << pairs_count_before - next_result.task_id() << ", ";
      estats_out << pairs_count_before << ", ";
      estats_out << next_result.s_.graph_size_ << ", ";
      estats_out << next_result.s_.edges_before_reweight_ << ", ";
      estats_out << next_result.s_.edges_after_reweight_ << ", ";
      estats_out << next_result.s_.count_after_harvest_ << ", ";
      estats_out << next_result.s_.count_after_normalize_ << ", ";
      estats_out << pairs_added << ", ";
      estats_out << next_result.u().size() << ", ";
      estats_out << next_result.v().size() << ", ";
      estats_out << next_result.s_.worker_id_ << ", ";
      estats_out << total_time.folding().last() << ", ";
      estats_out << total_time.folding().average() << ", ";
      estats_out << total_time.reweight().last() << ", ";
      estats_out << total_time.reweight().average() << ", ";
      estats_out << total_time.harvest().last() << ", ";
      estats_out << total_time.harvest().average() << ", ";
      estats_out << total_time.normalize().last() << ", ";
      estats_out << total_time.normalize().average() << "\n";
    }

    if (proc_words.is_open()) {
      proc_words << next_result.task_id() << ", ";
      PrintWord(next_result.u(), &proc_words);
      proc_words << ", ";
      PrintWord(next_result.v(), &proc_words);
      proc_words << "\n";
    }

    out->flush();
    estats_out.flush();
    proc_words.flush();
    unproc_words.flush();

    future_results.pop_front();
  }
  overall_time_stamp.Click();
  overall_time_stamp.Reset();
  auto average_other_time = (overall_time.total_time() - total_time.total()) / total_time.folding().iterations();
  *out << "total time: " << overall_time.total_time() << std::endl;
  *out << "total rest(av): " << average_other_time << std::endl;
  *out << "total iterations: " << total_time.folding().iterations() << std::endl;

  return 0;
}
