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

   g++ -O3 -std=c++11 -Wall -Wextra -pedantic -o kauffman kauffman.cc

On linux you might have to do:

   g++ -O3 -pthread -std=c++11 -Wall -Wextra -pedantic -o kauffman kauffman.cc

*******************************************************************************/

#include <math.h>

#include <atomic>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <mutex>
#include <stack>
#include <thread>
#include <vector>

#include "Dyck/dyck.h"
#include "timer.h"

static const size_t max_nr_threads = std::thread::hardware_concurrency() - 2;

typedef uint_fast8_t          letter_t;
typedef std::vector<letter_t> dyck_word_t;
typedef size_t                dyck_index_t;

static std::mutex               mtx;
static bool                     verbose;

static std::vector<dyck_word_t>       DYCK_WORDS;
static std::vector<std::vector<bool>> DYCK_OUTER_BOOL;

static size_t const catalan_numbers[] = {1,
                                         1,
                                         2,
                                         5,
                                         14,
                                         42,
                                         132,
                                         429,
                                         1430,
                                         4862,
                                         16796,
                                         58786,
                                         208012,
                                         742900,
                                         2674440,
                                         9694845,
                                         35357670,
                                         129644790,
                                         477638700,
                                         1767263190,
                                         6564120420,
                                         24466267020,
                                         91482563640,
                                         343059613650,
                                         1289904147324,
                                         4861946401452,
                                         18367353072152,
                                         69533550916004,
                                         263747951750360,
                                         1002242216651368,
                                         3814986502092304};

void print_help_and_exit(char* name) {
  std::cout << "usage: " << name << " [-h] [-v] n" << std::endl;
  exit(0);
}

template <typename T> void print_vector(std::vector<T> const& v) {
  for (auto it = v.begin(); it != v.end(); ++it) {
    std::cout << *it << ", ";
  }
  std::cout << std::endl;
}

void print_mem_usage(Timer& timer) {
  timer.print();
  std::cout << std::endl;

  double mem = 0;

  size_t len = DYCK_WORDS[0].size();

  mem += (DYCK_WORDS.size() * sizeof(letter_t) * len);
  mem += (DYCK_WORDS.size() * sizeof(bool) * len);

  std::string suf;
  if (mem > 1073741824) {  // 1024 ^ 3
    mem = mem / 1073741824;
    suf = " GB";
  } else if (mem > 1048576) {  // 1024 ^ 2
    mem = mem / 1048576;
    suf = " MB";
  } else if (mem > 1024) {
    mem = mem / 1024;
    suf = " KB";
  } else {
    suf = " bytes";
  }
  std::cout << "Dyck words use ~ " << mem << suf << std::endl;
  std::cout << "Using " << max_nr_threads << " / "
            << std::thread::hardware_concurrency() << " threads" << std::endl;
}

void count_even(size_t       deg,
                size_t       thread_id,
                dyck_index_t nr_dyck_words,
                size_t       begin,
                size_t       end,
                size_t&      nr_idempotents) {
  Timer timer;
  if (verbose) {
    timer.start();
  }

  std::vector<bool> seen(deg, false);

  for (dyck_index_t i = begin; i < end; i++) {
    for (dyck_index_t j = i + 1; j < nr_dyck_words; j++) {
      std::fill(seen.begin(), seen.end(), false);
      size_t cnt = 1, pos = 0;
      while (pos < deg) {
        size_t nr_i = 0, nr_j = 0;

        if (DYCK_OUTER_BOOL[i][pos]) {
          nr_i++;
        }
        if (DYCK_OUTER_BOOL[j][pos]) {
          nr_j++;
        }
        seen[pos]                = true;
        seen[DYCK_WORDS[j][pos]] = true;
        pos                      = DYCK_WORDS[i][DYCK_WORDS[j][pos]];

        while (!seen[pos]) {
          seen[pos]                = true;
          seen[DYCK_WORDS[j][pos]] = true;
          if (DYCK_OUTER_BOOL[j][pos]) {
            nr_j++;
          } else if (DYCK_OUTER_BOOL[i][pos]) {
            nr_i++;
            pos = DYCK_WORDS[i][DYCK_WORDS[j][pos]];
            break;
          }
          pos = DYCK_WORDS[i][DYCK_WORDS[j][pos]];
        }
        while (!seen[pos]) {
          if (DYCK_OUTER_BOOL[i][pos]) {
            nr_i++;
          }
          seen[pos]                = true;
          seen[DYCK_WORDS[j][pos]] = true;
          pos                      = DYCK_WORDS[i][DYCK_WORDS[j][pos]];
        }
        if (nr_i == 0 || nr_j == 0) {
          cnt = 0;
          break;
        }
        cnt *= (nr_i * nr_j);
        while (seen[pos]) pos++;
      }
      if (cnt != 0) {
        nr_idempotents += (2 * cnt);
      }
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

void count_odd(size_t       dyck_word_length,
               size_t       thread_id,
               dyck_index_t nr_dyck_words,
               size_t       begin,
               size_t       end,
               size_t&      nr_idempotents) {
  Timer timer;
  if (verbose) {
    timer.start();
  }
  assert(dyck_word_length == DYCK_WORDS[0].size());
  std::vector<bool> seen(dyck_word_length, false);

  for (dyck_index_t i = begin; i < end; i++) {
    for (dyck_index_t j = i + 1; j < nr_dyck_words; j++) {
      // check if the extra point is in a path that contains all of the
      // elements
      std::fill(seen.begin(), seen.end(), false);
      size_t pos = dyck_word_length - 1, nr_seen = 0;
      size_t cutoff = dyck_word_length;

      while (!seen[pos]) {
        nr_seen += 2;
        seen[pos] = true;
        pos       = DYCK_WORDS[j][pos];
        seen[pos] = true;
        cutoff    = (pos < cutoff ? pos : cutoff);
        if (DYCK_OUTER_BOOL[i][pos]) {
          pos = DYCK_WORDS[i][pos];
          break;
        }
        pos = DYCK_WORDS[i][pos];
      }
      while (!seen[pos]) {
        nr_seen += 2;
        seen[pos] = true;
        pos       = DYCK_WORDS[j][pos];
        seen[pos] = true;
        pos       = DYCK_WORDS[i][pos];
      }
      if (nr_seen == dyck_word_length) {
        nr_idempotents += 2;
        continue;
      }

      size_t cnt = 1;
      pos        = 0;

      while (pos < cutoff) {
        size_t nr_i = 0, nr_j = 0;

        if (DYCK_OUTER_BOOL[i][pos]) {
          nr_i++;
        }
        if (DYCK_OUTER_BOOL[j][pos]) {
          nr_j++;
        }
        seen[pos]                = true;
        seen[DYCK_WORDS[j][pos]] = true;
        pos                      = DYCK_WORDS[i][DYCK_WORDS[j][pos]];

        while (!seen[pos]) {
          seen[pos]                = true;
          seen[DYCK_WORDS[j][pos]] = true;
          if (DYCK_OUTER_BOOL[j][pos]) {
            nr_j++;
          } else if (DYCK_OUTER_BOOL[i][pos]) {
            nr_i++;
            pos = DYCK_WORDS[i][DYCK_WORDS[j][pos]];
            break;
          }
          pos = DYCK_WORDS[i][DYCK_WORDS[j][pos]];
        }
        while (!seen[pos]) {
          if (DYCK_OUTER_BOOL[i][pos]) {
            nr_i++;
          }
          seen[pos]                = true;
          seen[DYCK_WORDS[j][pos]] = true;
          pos                      = DYCK_WORDS[i][DYCK_WORDS[j][pos]];
        }
        if (nr_i == 0 || nr_j == 0) {
          cnt = 0;
          break;
        }
        cnt *= (nr_i * nr_j);
        while (seen[pos] && pos <= cutoff) pos++;
      }

      // pos > cutoff
      if (cnt != 0) {
        for (; pos < dyck_word_length; pos++) {
          if (!seen[pos]) {
            cnt = 0;
            break;
          }
        }
        if (cnt != 0) {
          nr_idempotents += (2 * cnt);
        }
      }
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
  // Not very robust parsing!
  long deg = 0;
  verbose  = false;

  while (--argc) {
    if (argv[argc][0] == '-') {

      const char* p = argv[argc];

      while (*++p) {

        switch (*p) {

          case 'v':
            verbose = true;
            break;

          case 'h':
            print_help_and_exit(argv[0]);
        }
      }
    } else {
      char** end = nullptr;
      deg        = strtol(argv[argc], end, 0);
      if (deg <= 0 || deg > 40) {
        std::cerr << argv[0] << ": invalid argument! " << std::endl
                  << argv[0] << ": must be an integer in [1, 30]" << std::endl;
        exit(-1);
      }
    }
  }

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

  DYCK_WORDS.reserve(nr_dyck_words);
  DYCK_OUTER_BOOL.reserve(nr_dyck_words);

  dyck::integer            mask;
  std::stack<dyck_index_t> stack;
  ;
  dyck::integer w = dyck::minimum(n);

  for (dyck_index_t i = 0; i < nr_dyck_words; i++, w = dyck::next(w)) {
    mask = static_cast<dyck::integer>(1) << (2 * n - 1);
    DYCK_WORDS.push_back(dyck_word_t());
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
    DYCK_OUTER_BOOL.push_back(std::vector<bool>(2 * n, false));
    for (dyck_index_t j = 0; j < 2 * n; j = DYCK_WORDS[i][j], j++) {
      DYCK_OUTER_BOOL[i][j] = true;
    }
  }

  if (verbose) {
    print_mem_usage(timer);
  }
  size_t out = 1;
  if (nr_dyck_words < 400) {
    if ((deg / 2) * 2 == deg) {  // deg is even
      count_even(2 * n, 0, nr_dyck_words, 0, nr_dyck_words, out);
    } else {
      count_odd(2 * n, 0, nr_dyck_words, 0, nr_dyck_words, out);
    }
  } else {

    size_t av_load =
        (nr_dyck_words * (nr_dyck_words - 1)) / (2 * max_nr_threads);
    size_t              thread_id   = 0;
    size_t              thread_load = 0;
    std::vector<size_t> begin(max_nr_threads, 0);
    std::vector<size_t> end(max_nr_threads, nr_dyck_words);

    for (dyck_index_t i = 0; i < nr_dyck_words; i++) {
      thread_load += nr_dyck_words - i - 1;
      if (thread_load >= av_load && thread_id != max_nr_threads - 1) {
        end[thread_id] = i + 1;
        thread_id++;
        begin[thread_id] = i + 1;
        thread_load      = 0;
      }
    }
    size_t                   nr_threads = thread_id + 1;
    std::vector<size_t>      nr_idempotents(nr_threads, 0);
    std::vector<std::thread> threads;

    if ((deg / 2) * 2 == deg) {  // deg is even
      for (size_t i = 0; i < nr_threads; i++) {
        threads.push_back(std::thread(count_even,
                                      2 * n,
                                      i,
                                      nr_dyck_words,
                                      begin[i],
                                      end[i],
                                      std::ref(nr_idempotents[i])));
      }
    } else {
      for (size_t i = 0; i < nr_threads; i++) {
        threads.push_back(std::thread(count_odd,
                                      2 * n,
                                      i,
                                      nr_dyck_words,
                                      begin[i],
                                      end[i],
                                      std::ref(nr_idempotents[i])));
      }
    }

    for (size_t i = 0; i < nr_threads; i++) {
      threads[i].join();
      out += nr_idempotents[i];
    }
  }

  if (verbose) {
    std::cout << "Total elapsed time = ";
    timer.print();
    std::cout << std::endl;
  }
  std::cout << out << std::endl;
  exit(0);
}
