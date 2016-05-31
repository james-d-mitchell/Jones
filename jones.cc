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
#include <stack>
#include <thread>
#include <mutex>
#include <vector>

#include "Dyck/dyck.h"
#include "timer.h"

typedef std::vector<u_int8_t> dyck_word_t;
typedef size_t                dyck_index_t;

static const size_t nr_threads = std::thread::hardware_concurrency() - 2;
static std::mutex   mtx;
static bool         verbose;

static std::vector<dyck_word_t>       dyck_words;
static std::vector<dyck_word_t>       dyck_outer;
static std::vector<std::vector<bool>> dyck_outer_memb;

static size_t const catalan_numbers[] = {1, 1, 2, 5, 14, 42, 132, 429, 1430,
  4862, 16796, 58786, 208012, 742900, 2674440, 9694845, 35357670, 129644790,
  477638700, 1767263190, 6564120420, 24466267020, 91482563640, 343059613650,
  1289904147324, 4861946401452, 18367353072152, 69533550916004,
  263747951750360, 1002242216651368, 3814986502092304};

void print_help_and_exit(char* name) {
  std::cout << "usage: " << name << " [-h] [-v] n" << std::endl;
  exit(0);
}

void count_even (size_t                           thread_id,
                 dyck_index_t                     nr_dyck_words,
                 std::vector<dyck_index_t> const& unprocessed,
                 size_t&                          nr_idempotents) {
  Timer timer;
  if (verbose) {
    timer.start();
  }
  for (dyck_index_t i: unprocessed) {
    nr_idempotents += pow(2, dyck_outer[i].size());
    for (dyck_index_t j = i + 1; j < nr_dyck_words; j++) {
      size_t max = 0, cnt = 1;
      dyck_word_t::iterator it = dyck_outer[j].begin();
      do {
        while (*it < max) {
          it++;
        }
        size_t pos = *it;
        size_t nr_i = 0, nr_j = 1;

        max = dyck_words[j][pos];

        if (dyck_outer_memb[i][pos]) {
          nr_i++;
        }
        pos = dyck_words[i][dyck_words[j][pos]];

        while (*it != pos) {
          if (dyck_outer_memb[j][pos]) {
            nr_j++;
            if (dyck_words[j][pos] > max) {
              max = dyck_words[j][pos];
            }
          } else if (dyck_outer_memb[i][pos]) {
            nr_i++;
            pos = dyck_words[i][dyck_words[j][pos]];
            break;
          }
          pos = dyck_words[i][dyck_words[j][pos]];
        }
        while (*it != pos) {
          if (dyck_outer_memb[i][pos]) {
            nr_i++;
          }
          pos = dyck_words[i][dyck_words[j][pos]];
        }
        cnt *= (nr_i * nr_j + 1);
      } while (max < dyck_outer[j].back());
      nr_idempotents += (2 * cnt);
    }
  }
  if (verbose) {
    mtx.lock();
    std::cout << "Thread " << thread_id << " is finished, ";
    timer.stop();
    mtx.unlock();
  }
}

void count_odd (size_t                            thread_id,
                dyck_index_t                      nr_dyck_words,
                std::vector<dyck_index_t> const&  unprocessed,
                size_t&                           nr_idempotents) {
  Timer timer;
  if (verbose) {
    timer.start();
  }
  for (dyck_index_t i: unprocessed) {
    size_t n = dyck_words[i][dyck_words[0].size() - 1];

    nr_idempotents += pow(2, dyck_outer[i].size() - 1);
    for (dyck_index_t j = i + 1; j < nr_dyck_words; j++) {
      size_t max = 0, cnt = 1, pos;
      dyck_word_t::iterator it = dyck_outer[j].begin();
      do {
        while (*it < max) {
          it++;
        }
        size_t nr_i = 0, nr_j = 1;
        max = dyck_words[j][*it];
        if (dyck_outer_memb[i][*it]) {
          nr_i++;
        }
        pos = dyck_words[i][dyck_words[j][*it]];

        while (pos != *it && pos != n) {
          if (dyck_outer_memb[j][pos]) {
            nr_j++;
            if (dyck_words[j][pos] > max) {
              max = dyck_words[j][pos];
            }
          } else if (dyck_outer_memb[i][pos]) {
            nr_i++;
            pos = dyck_words[i][dyck_words[j][pos]];
            break;
          }
          pos = dyck_words[i][dyck_words[j][pos]];
        }
        while (pos != *it && pos != n) {
          if (dyck_outer_memb[i][pos]) {
            nr_i++;
          }
          pos = dyck_words[i][dyck_words[j][pos]];
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
    std::cout << "Thread " << thread_id << " is finished, ";
    timer.stop();
    mtx.unlock();
  }
}

int main(int argc, char* argv[]) {
  // Not very robust parsing!
  long deg = 0;
  verbose = false;

  while (--argc) {
    if (argv[argc][0] == '-') {

      const char* p = argv[argc];

      while (*++p) {

        switch (*p) {

          case 'v' :
            verbose = true;
            break;

          case 'h' :
            print_help_and_exit(argv[0]);
        }
      }
    } else {
      char** end = nullptr;
      deg = strtol(argv[argc], end, 0);
      if (deg <= 0 || deg > 40) {
        std::cerr <<
          argv[0] << ": invalid argument! " << std::endl <<
          argv[0] << ": must be an integer in [1, 30]" << std::endl;
        exit(-1);
      }
    }
  }

  if (deg == 0) {
    print_help_and_exit(argv[0]);
  }

  dyck_index_t n;
  if ((deg / 2) * 2 == deg) { // deg is even
    n = deg / 2; // input to dyck, half the length of the returned words
  } else {
    n = (deg + 1) / 2;
  }
  size_t const nr_dyck_words = catalan_numbers[n];

  Timer timer;
  if (verbose) {
    std::cout << "Number of Dyck words is " << nr_dyck_words << std::endl;
    std::cout << "Processing Dyck words, ";
    timer.start();
  }

  dyck_words = std::vector<dyck_word_t>();
  dyck_words.reserve(nr_dyck_words);

  dyck_outer = std::vector<dyck_word_t>();
  dyck_outer.reserve(nr_dyck_words);

  dyck_outer_memb = std::vector<std::vector<bool>>();
  dyck_outer_memb.reserve(nr_dyck_words);

  dyck::integer mask;
  std::stack<dyck_index_t> stack = std::stack<dyck_index_t>();
  dyck::integer w = dyck::minimum(n);

  for (dyck_index_t i = 0; i < nr_dyck_words; i++, w = dyck::next(w)) {
    mask = static_cast<dyck::integer>(1) << (2 * n - 1);
    dyck_words.push_back(dyck_word_t());
    dyck_words[i].resize(2 * n);

    for (dyck_index_t j = 0; j < 2 * n; j++, mask >>= 1) {
      if (mask & w) { // opening bracket
        stack.push(j);
      } else {
        dyck_words[i][j] = stack.top();
        dyck_words[i][stack.top()] = j;
        stack.pop();
      }
    }
    dyck_outer.push_back(dyck_word_t());
    dyck_outer_memb.push_back(std::vector<bool>(2 * n, false));
    for (dyck_index_t j = 0; j < 2 * n; j = dyck_words[i][j], j++) {
      dyck_outer[i].push_back(j);
      dyck_outer_memb[i][j] = true;
    }
  }

  if (verbose) {
    timer.print();
    double mem = 0;

    for (auto x: dyck_outer) {
      mem += (x.size() * sizeof(u_int8_t));
    }
    mem += (nr_dyck_words * sizeof(u_int8_t) * 2 * n);
    mem += (nr_dyck_words * sizeof(bool) * 2 * n);

    std::string suf;
    if (mem > 1073741824) { // 1024 ^ 3
      mem = mem / 1073741824;
      suf = " GB";
    } else if (mem > 1048576) { // 1024 ^ 2
      mem = mem / 1048576;
      suf = " MB";
    } else if (mem > 1024) {
      mem = mem / 1024;
      suf = " KB";
    } else {
      suf = " bytes";
    }

    std::cout << "Dyck words use ~ " << mem << suf << std::endl;

    std::cout << "Using " << nr_threads << " / " <<
      std::thread::hardware_concurrency() << " threads" << std::endl;
  }

  size_t av_load     = (nr_dyck_words * (nr_dyck_words - 1))
                        / (2 * nr_threads);
  size_t thread_id   = 0;
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
  std::vector<std::thread> threads = std::vector<std::thread>();

  if ((deg / 2) * 2 == deg) { // deg is even
    for (size_t i = 0; i < nr_threads; i++) {
      threads.push_back(std::thread(count_even,
                                    i,
                                    nr_dyck_words,
                                    std::ref(unprocessed[i]),
                                    std::ref(nr_idempotents[i])));
    }
  } else {
    for (size_t i = 0; i < nr_threads; i++) {
      threads.push_back(std::thread(count_odd,
                                    i,
                                    nr_dyck_words,
                                    std::ref(unprocessed[i]),
                                    std::ref(nr_idempotents[i])));
    }
  }

  size_t out = 0;
  for (size_t i = 0; i < nr_threads; i++) {
    threads[i].join();
    out += nr_idempotents[i];
  }

  if (verbose) {
    std::cout << "Total ";
    timer.stop();
  }
  std::cout << out << std::endl;
  exit(0);
}
