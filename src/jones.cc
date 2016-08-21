/*******************************************************************************

 Copyright (C) 2016 James D. Mitchell

 This work is licensed under a Creative Commons Attribution-ShareAlike 4.0
 International License. See
 http://creativecommons.org/licenses/by-sa/4.0/

 For a discussion about this code and latest version:

 This requires the code from the git repository:

   https://github.com/cassioneri/Dyck

 to be put inside a subdirectory (of the directory where this will be compiled)
 named Dyck.

 Compile with:

   g++ -O3 -std=c++11 -Wall -Wextra -pedantic -o jones jones.cc

On linux you might have to do:

   g++ -O3 -pthread -std=c++11 -Wall -Wextra -pedantic -o jones jones.cc

*******************************************************************************/

#include <math.h>

#include <cstdlib>
#include <functional>
#include <iostream>
#include <mutex>
#include <numeric>
#include <stack>
#include <thread>
#include <unordered_set>
#include <vector>

#include "base.h"

// User definable globals
static const size_t nr_threads = std::thread::hardware_concurrency() - 2;

// Simple class to contain the required information about a Dyck word

class Dyck {
 public:
  explicit Dyck(dyck::integer w, size_t n) : word(), outer(), lookup() {
    word.resize(2 * n, 0);
    lookup.resize(2 * n, false);

    dyck::integer mask = static_cast<dyck::integer>(1) << (2 * n - 1);

    for (dyck_index_t j = 0; j < 2 * n; j++, mask >>= 1) {
      if (mask & w) {  // opening bracket
        stack.push(j);
      } else {
        word[j]           = stack.top();
        word[stack.top()] = j;
        if (stack.size() == 1) {  // stack.top() is an outer bracket
          outer.push_back(stack.top());
          lookup[stack.top()] = true;
        }
        stack.pop();
      }
    }
  }

  dyck_vec_t        word;
  dyck_vec_t        outer;
  std::vector<bool> lookup;

  static std::stack<dyck_index_t> stack;

  inline letter_t const& operator()(letter_t const& i) const {
    return word[i];
  }
};

std::stack<dyck_index_t> Dyck::stack;

// Globals
static std::mutex mtx;
static bool       verbose;

// Dycks
static std::vector<Dyck> PALIN;
static std::vector<Dyck> NONPALIN;
static std::vector<Dyck> NONPALIN_R;

static std::vector<dyck_vec_t>        DYCK_WORDS;
static std::vector<dyck_vec_t>        DYCK_OUTER;
static std::vector<std::vector<bool>> DYCK_BOOL;

// Utility functions
void print_mem_usage_even(size_t deg) {
  double mem = (deg * catalan_numbers[deg / 2] * sizeof(letter_t));
  mem += (deg * catalan_numbers[deg / 2] * sizeof(bool));

  for (auto const& x : PALIN) {
    mem += (x.outer.size() * sizeof(letter_t));
  }
  for (auto const& x : NONPALIN) {
    mem += (x.outer.size() * sizeof(letter_t));
  }
  for (auto const& x : NONPALIN_R) {
    mem += (x.outer.size() * sizeof(letter_t));
  }

  std::cout << "Dyck words use ~ " << string_mem(mem) << std::endl;
  std::cout << "Using " << nr_threads << " / "
            << std::thread::hardware_concurrency() << " threads" << std::endl;
}

void print_mem_usage_odd(size_t n) {
  double mem = ((2 * n + 1) * catalan_numbers[n] * sizeof(letter_t));
  mem += ((2 * n + 1) * catalan_numbers[n] * sizeof(bool));

  for (auto const& x : DYCK_OUTER) {
    mem += (x.size() * sizeof(letter_t));
  }

  std::cout << "Dyck words use ~ " << string_mem(mem) << std::endl;
  std::cout << "Using " << nr_threads << " / "
            << std::thread::hardware_concurrency() << " threads" << std::endl;
}

void print_thread_finished(size_t thread_id, Timer& timer) {
  if (verbose) {
    mtx.lock();
    std::cout << "Thread " << thread_id << " is finished, elapsed time = ";
    timer.print();
    std::cout << std::endl;
    mtx.unlock();
  }
}

// The main event from lower to higher level

inline void count_cycle(size_t&     nr_idempotents,
                        size_t      multiplier,
                        Dyck const& u,
                        Dyck const& l) {
  size_t                     max = 0, cnt = 1;
  dyck_vec_t::const_iterator it = l.outer.cbegin();
  do {
    while (*it < max) it++;
    size_t pos  = *it;
    size_t nr_u = 0, nr_l = 1;

    max = l(pos);

    if (u.lookup[pos]) nr_u++;

    pos = u(l(pos));

    while (*it != pos) {
      if (l.lookup[pos]) {
        nr_l++;
        max = l(pos);
      } else if (u.lookup[pos]) {
        nr_u++;
        pos = u(l(pos));
        break;
      }
      pos = u(l(pos));
    }
    while (*it != pos) {
      if (u.lookup[pos]) nr_u++;
      pos = u(l(pos));
    }
    cnt *= (nr_u * nr_l + 1);
  } while (max < l.outer.back());
  nr_idempotents += (multiplier * cnt);
}

void count_even_tri_thread(size_t                           thread_id,
                           std::vector<dyck_index_t> const& index,
                           std::vector<Dyck> const&         dycks1,
                           std::vector<Dyck> const&         dycks2,
                           size_t&                          nr_idempotents,
                           size_t                           multiplier) {
  // dycks2 is only here to make the signature of the function correct
  (void) dycks2;
  Timer timer;
  if (verbose) {
    timer.start();
  }
  size_t nr = dycks1.size();
  for (auto const& i : index) {
    for (size_t j = i + 1; j < nr; j++) {
      count_cycle(nr_idempotents, multiplier, dycks1[i], dycks1[j]);
    }
  }
  print_thread_finished(thread_id, timer);
}

void count_even_rect_thread(size_t                           thread_id,
                            std::vector<dyck_index_t> const& index,
                            std::vector<Dyck> const&         dycks1,
                            std::vector<Dyck> const&         dycks2,
                            size_t&                          nr_idempotents,
                            size_t                           multiplier) {
  Timer timer;
  if (verbose) {
    timer.start();
  }

  for (auto const& i : index) {
    for (size_t j = 0; j < dycks2.size(); j++) {
      count_cycle(nr_idempotents, multiplier, dycks1[i], dycks2[j]);
    }
  }
  print_thread_finished(thread_id, timer);
}

void count_even_reverse_thread(size_t                           thread_id,
                               std::vector<dyck_index_t> const& index,
                               std::vector<Dyck> const&         dycks1,
                               std::vector<Dyck> const&         dycks2,
                               size_t&                          nr_idempotents,
                               size_t                           multiplier) {
  Timer timer;
  if (verbose) {
    timer.start();
  }
  assert(multiplier == 4);
  size_t nr_dyck2 = dycks2.size();
  for (auto const& i : index) {
    count_cycle(nr_idempotents, 2, dycks1[i], dycks2[i]);
    for (size_t j = i + 1; j < nr_dyck2; j++) {
      count_cycle(nr_idempotents, 4, dycks1[i], dycks2[j]);
    }
  }
  print_thread_finished(thread_id, timer);
}

void distribute_to_threads(std::vector<Dyck> const&      dycks1,
                           std::vector<Dyck> const&      dycks2,
                           std::vector<size_t>&          nr_idempotents,
                           size_t                        multiplier,
                           size_t const                  av_load,
                           std::function<size_t(size_t)> cost,
                           std::function<void(size_t,
                                              std::vector<dyck_index_t> const&,
                                              std::vector<Dyck> const&,
                                              std::vector<Dyck> const&,
                                              size_t&,
                                              size_t)> thread_func) {
  std::vector<std::thread>               threads;
  std::vector<std::vector<dyck_index_t>> index(nr_threads,
                                               std::vector<dyck_index_t>());

  size_t const nr_dycks    = dycks1.size();
  size_t       thread_id   = 0;
  size_t       thread_load = 0;

  for (size_t i = 0; i < nr_dycks; i++) {
    index[thread_id].push_back(i);
    thread_load += cost(i);
    if (thread_load >= av_load && thread_id != nr_threads - 1) {
      if (verbose) {
        std::cout << "Thread " << thread_id << " has load " << thread_load
                  << std::endl;
      }
      thread_id++;
      thread_load = 0;
    }
  }
  if (verbose) {
    std::cout << "Thread " << thread_id << " has load " << thread_load
              << std::endl;
  }

  for (size_t i = 0; i < nr_threads; i++) {
    threads.push_back(std::thread(thread_func,
                                  i,
                                  std::ref(index[i]),
                                  std::ref(dycks1),
                                  std::ref(dycks2),
                                  std::ref(nr_idempotents[i]),
                                  multiplier));
  }
  for (size_t i = 0; i < nr_threads; i++) {
    threads[i].join();
  }
}

void count_even_tri(std::vector<Dyck> const& dycks,
                    std::vector<size_t>&     nr_idempotents,
                    size_t const             multiplier) {
  size_t const nr      = dycks.size();
  size_t const av_load = (nr * (nr - 1)) / (2 * nr_threads);

  distribute_to_threads(dycks,
                        dycks,
                        nr_idempotents,
                        multiplier,
                        av_load,
                        [nr](size_t i) { return nr - i - 1; },
                        count_even_tri_thread);
}

void count_even_rect(std::vector<Dyck> const& dycks1,
                     std::vector<Dyck> const& dycks2,
                     std::vector<size_t>&     nr_idempotents) {
  size_t const nr_dycks1 = dycks1.size();
  size_t const nr_dycks2 = dycks2.size();
  size_t const av_load   = (nr_dycks1 * nr_dycks2) / nr_threads;

  distribute_to_threads(dycks1,
                        dycks2,
                        nr_idempotents,
                        4,
                        av_load,
                        [nr_dycks2](size_t i) {
                          (void) i;
                          return nr_dycks2;
                        },
                        count_even_rect_thread);
}

void count_even_reverse(std::vector<Dyck> const& dycks1,
                        std::vector<Dyck> const& dycks2,
                        std::vector<size_t>&     nr_idempotents) {
  size_t const nr_dycks1 = dycks1.size();
  size_t const nr_dycks2 = dycks2.size();
  size_t const av_load   = (nr_dycks1 * (nr_dycks2 + 1)) / (2 * nr_threads);

  distribute_to_threads(dycks1,
                        dycks2,
                        nr_idempotents,
                        4,
                        av_load,
                        [nr_dycks2](size_t i) { return nr_dycks2 - i; },
                        count_even_reverse_thread);
}

// The code for the odd case is simpler but involves doing ~2 times more
// comparisons

void count_odd(size_t                           thread_id,
               dyck_index_t                     nr_dyck_words,
               std::vector<dyck_index_t> const& unprocessed,
               size_t&                          nr_idempotents) {
  Timer timer;
  if (verbose) {
    timer.start();
  }
  for (dyck_index_t i : unprocessed) {
    size_t n = DYCK_WORDS[i][DYCK_WORDS[0].size() - 1];

    nr_idempotents += pow(2, DYCK_OUTER[i].size() - 1);
    for (dyck_index_t j = i + 1; j < nr_dyck_words; j++) {
      size_t               max = 0, cnt = 1, pos;
      dyck_vec_t::iterator it = DYCK_OUTER[j].begin();
      do {
        while (*it < max) {
          it++;
        }
        size_t nr_i = 0, nr_j = 1;
        max = DYCK_WORDS[j][*it];
        if (DYCK_BOOL[i][*it]) {
          nr_i++;
        }
        pos = DYCK_WORDS[i][DYCK_WORDS[j][*it]];

        while (pos != *it && pos != n) {
          if (DYCK_BOOL[j][pos]) {
            nr_j++;
            max = DYCK_WORDS[j][pos];
          } else if (DYCK_BOOL[i][pos]) {
            nr_i++;
            pos = DYCK_WORDS[i][DYCK_WORDS[j][pos]];
            break;
          }
          pos = DYCK_WORDS[i][DYCK_WORDS[j][pos]];
        }
        while (pos != *it && pos != n) {
          if (DYCK_BOOL[i][pos]) {
            nr_i++;
          }
          pos = DYCK_WORDS[i][DYCK_WORDS[j][pos]];
        }
        if (pos != n) {
          cnt *= (nr_j * nr_i + 1);
        }
      } while (pos != n);
      nr_idempotents += (2 * cnt);
    }
  }
  if (verbose) {
    mtx.lock();
    std::cout << "Thread " << thread_id << " is finished, elapsed time = ";
    timer.print();
    std::cout << std::endl;
    mtx.unlock();
  }
}
int main(int argc, char* argv[]) {
  verbose    = false;
  size_t deg = 0;

  parse_args(argc, argv, verbose, deg);

  if (deg == 0) {
    print_help_and_exit(argv[0]);
  }

  dyck_index_t n;
  if ((deg / 2) * 2 == deg) {  // deg is even
    n = deg / 2;  // input to dyck, half the length of the returned words
  } else {
    n = (deg + 1) / 2;
  }
  size_t const nr_dyck_words = catalan_numbers[n];

  Timer timer;
  if (verbose) {
    std::cout << "Number of Dyck words is " << nr_dyck_words << std::endl;
    std::cout << "Processing Dyck words, elapsed time = ";
    timer.start();
  }

  if ((deg / 2) * 2 == deg) {
    // Number of idempotents arising from (w, w):
    size_t palin    = 0;  // where w is a palindromic Dyck word
    size_t nonpalin = 0;  // where w is a non-palindromic Dyck word

    {
      dyck::integer                     w = dyck::minimum(n);
      std::unordered_set<dyck::integer> reversed;

      for (dyck_index_t i = 0; i < nr_dyck_words; i++, w = dyck::next(w)) {
        dyck::integer ww = reverse(w, 2 * n);
        Dyck          next(w, n);

        if (ww == w) {
          palin += pow(2, next.outer.size());
          PALIN.push_back(next);
        } else if (reversed.find(ww) == reversed.end()) {
          nonpalin += pow(2, next.outer.size() + 1);
          reversed.insert(w);
          NONPALIN.push_back(next);
          NONPALIN_R.push_back(Dyck(ww, n));
        }
      }
    }
    std::vector<size_t> nr_idempotents(nr_threads, 0);
    size_t              last = 0;

    if (verbose) {
      std::cout << timer.string() << std::endl;
      print_mem_usage_even(2 * n);
      std::cout << "Number of palindromic Dyck words is " << PALIN.size()
                << std::endl;
      std::cout << "Number of non-palindromic Dyck words is " << NONPALIN.size()
                << std::endl;
    }
    assert(NONPALIN.size() == NONPALIN_R.size());

    count_even_tri(PALIN, nr_idempotents, 2);
    if (verbose) {
      last = std::accumulate(nr_idempotents.begin(), nr_idempotents.end(), 0);
      std::cout << "From comparison of palindromic and palindromic: "
                << last + palin << std::endl;
    }

    count_even_tri(NONPALIN, nr_idempotents, 4);
    if (verbose) {
      size_t next = std::accumulate(
          nr_idempotents.begin(), nr_idempotents.end(), (size_t) 0);
      std::cout << "From comparison of non-palindromic and non-palindromic: "
                << next + nonpalin - last << std::endl;
      last = next;
    }

    count_even_rect(PALIN, NONPALIN, nr_idempotents);
    if (verbose) {
      size_t next = std::accumulate(
          nr_idempotents.begin(), nr_idempotents.end(), (size_t) 0);
      std::cout << "From comparison of palindromic and non-palindromic: "
                << next - last << std::endl;
      last = next;
    }

    count_even_reverse(NONPALIN, NONPALIN_R, nr_idempotents);
    if (verbose) {
      size_t next = std::accumulate(
          nr_idempotents.begin(), nr_idempotents.end(), (size_t) 0);
      std::cout << "From comparison of non-palindromics and their reverses: "
                << next - last << std::endl;
      std::cout << "Total elapsed time = " << timer.string() << std::endl;
    }
    std::cout << std::accumulate(
                     nr_idempotents.begin(), nr_idempotents.end(), (size_t) 0)
                     + palin + nonpalin
              << std::endl;
  } else {
    DYCK_WORDS.reserve(nr_dyck_words);
    DYCK_OUTER.reserve(nr_dyck_words);
    DYCK_BOOL.reserve(nr_dyck_words);

    dyck::integer            mask;
    std::stack<dyck_index_t> stack;
    dyck::integer w = dyck::minimum(n);

    for (dyck_index_t i = 0; i < nr_dyck_words; i++, w = dyck::next(w)) {
      mask = static_cast<dyck::integer>(1) << (2 * n - 1);
      DYCK_WORDS.push_back(dyck_vec_t());
      DYCK_WORDS[i].resize(2 * n);

      for (dyck_index_t j = 0; j < 2 * n; j++, mask >>= 1) {
        if (mask & w) {  // opening bracket
          stack.push(j);
        } else {
          DYCK_WORDS[i][j]           = stack.top();
          DYCK_WORDS[i][stack.top()] = j;
          stack.pop();
        }
      }
      DYCK_OUTER.push_back(dyck_vec_t());
      DYCK_BOOL.push_back(std::vector<bool>(2 * n, false));
      for (dyck_index_t j = 0; j < 2 * n; j = DYCK_WORDS[i][j], j++) {
        DYCK_OUTER[i].push_back(j);
        DYCK_BOOL[i][j] = true;
      }
    }
    if (verbose) {
      timer.print();
      std::cout << std::endl;
      print_mem_usage_odd(n);
    }
    size_t av_load   = (nr_dyck_words * (nr_dyck_words - 1)) / (2 * nr_threads);
    size_t thread_id = 0;
    size_t thread_load = 0;

    std::vector<std::vector<dyck_index_t>> unprocessed;
    for (size_t i = 0; i < nr_threads; i++) {
      unprocessed.push_back(std::vector<dyck_index_t>());
    }

    for (dyck_index_t i = 0; i < nr_dyck_words; i++) {
      unprocessed[thread_id].push_back(i);
      thread_load += nr_dyck_words - i - 1;
      if (thread_load >= av_load && thread_id != nr_threads - 1) {
        thread_id++;
        thread_load = 0;
      }
    }

    std::vector<size_t>      nr_idempotents(nr_threads, 0);
    std::vector<std::thread> threads;

    for (size_t i = 0; i < nr_threads; i++) {
      threads.push_back(std::thread(count_odd,
                                    i,
                                    nr_dyck_words,
                                    std::ref(unprocessed[i]),
                                    std::ref(nr_idempotents[i])));
    }

    size_t out = 0;
    for (size_t i = 0; i < nr_threads; i++) {
      threads[i].join();
      out += nr_idempotents[i];
    }

    if (verbose) {
      std::cout << "Total elapsed time = ";
      timer.print();
      std::cout << std::endl;
    }
    std::cout << out << std::endl;
  }
  exit(0);
}
