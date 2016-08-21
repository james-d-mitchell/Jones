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

   g++ -O3 -std=c++11 -Wall -Wextra -pedantic -o motzkin motzkin.cc

On linux you might have to do:

   g++ -O3 -pthread -std=c++11 -Wall -Wextra -pedantic -o motzkin motzkin.cc

*******************************************************************************/

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <mutex>
#include <random>
#include <stack>
#include <thread>
#include <vector>

#include "Dyck/dyck.h"
#include "timer.h"

typedef uint_fast8_t          letter_t;
typedef std::vector<letter_t> motzkin_word_t;
typedef size_t                index_t;
typedef uint32_t              subset_t;
typedef uint32_t              dyck_word_t;

static const size_t nr_threads = std::thread::hardware_concurrency() - 2;
static std::mutex   mtx;
static bool         verbose;

static std::vector<motzkin_word_t>    MOTZKIN_WORDS;
static std::vector<motzkin_word_t>    MOTZKIN_OUTER;
static std::vector<std::vector<bool>> MOTZKIN_BOOL;
static std::vector<dyck_word_t>       DYCK_WORDS;
static std::vector<subset_t>          SUBSETS;

static size_t const nr_motzkin_words_weight_0[] = {0,
                                                   1,
                                                   2,
                                                   4,
                                                   9,
                                                   21,
                                                   51,
                                                   127,
                                                   323,
                                                   835,
                                                   2188,
                                                   5798,
                                                   15511,
                                                   41835,
                                                   113634,
                                                   310572,
                                                   853467,
                                                   2356779,
                                                   6536382,
                                                   18199284,
                                                   50852019,
                                                   142547559,
                                                   400763223,
                                                   1129760415,
                                                   3192727797,
                                                   9043402501,
                                                   25669818476,
                                                   73007772802,
                                                   208023278209,
                                                   593742784829,
                                                   1697385471211,
                                                   4859761676391,
                                                   13933569346707,
                                                   40002464776083,
                                                   114988706524270,
                                                   330931069469828,
                                                   953467954114363,
                                                   2750016719520991,
                                                   7939655757745265,
                                                   22944749046030949,
                                                   66368199913921497};

// weight 1
static size_t const nr_motzkin_words_weight_1[] = {0,
                                                   1,
                                                   2,
                                                   5,
                                                   12,
                                                   30,
                                                   76,
                                                   196,
                                                   512,
                                                   1353,
                                                   3610,
                                                   9713,
                                                   26324,
                                                   71799,
                                                   196938,
                                                   542895,
                                                   1503312,
                                                   4179603,
                                                   11662902,
                                                   32652735,
                                                   91695540,
                                                   258215664,
                                                   728997192,
                                                   2062967382,
                                                   5850674704,
                                                   16626415975,
                                                   47337954326,
                                                   135015505407,
                                                   385719506620,
                                                   1103642686382,
                                                   3162376205180,
                                                   9073807670316,
                                                   26068895429376,
                                                   74986241748187,
                                                   215942362945558,
                                                   622536884644535,
                                                   1796548765406628,
                                                   5189639038224274,
                                                   15005093288285684,
                                                   43423450867890548,
                                                   125769718187920320};

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

void init_motzkin(size_t                        nr_motzkin_words,
                  size_t                        motzkin_word_length,
                  size_t                        dyck_length_min,
                  size_t                        dyck_length_max,
                  size_t                        set_size,
                  std::function<size_t(size_t)> subset_size) {
  MOTZKIN_WORDS.clear();
  MOTZKIN_WORDS.reserve(nr_motzkin_words);
  MOTZKIN_OUTER.clear();
  MOTZKIN_OUTER.reserve(nr_motzkin_words);
  MOTZKIN_BOOL.clear();
  MOTZKIN_BOOL.reserve(nr_motzkin_words);

  index_t i = 0;

  for (size_t m = dyck_length_min; m <= dyck_length_max; m++) {
    DYCK_WORDS.clear();
    SUBSETS.clear();
    // Get Dyck words of length m
    dyck::integer const end = dyck::maximum(m);
    for (dyck::integer w = dyck::minimum(m); w <= end; w = dyck::next(w)) {
      DYCK_WORDS.push_back(w);
    }

    size_t K = subset_size(m);
    assert(K + 2 * m == motzkin_word_length);

    if (K == 0) {
      SUBSETS.push_back(0);
    } else {
      // Get subsets of {1 .. set_size} of size subset_size(m)
      int s = (1 << K) - 1;
      while (!(s & 1 << set_size)) {
        SUBSETS.push_back(s);
        int lo = s & ~(s - 1);   // lowest one bit
        int lz = (s + lo) & ~s;  // lowest zero bit above lo
        s |= lz;                 // add lz to the set
        s &= ~(lz - 1);          // reset bits below lz
        s |= (lz / lo / 2) - 1;  // put back right number of bits at end
      }
    }

    dyck::integer       mask_word;
    dyck::integer       mask_subset;
    std::stack<index_t> stack;

    for (dyck_word_t const& w : DYCK_WORDS) {
      for (subset_t const& s : SUBSETS) {
        mask_word   = static_cast<dyck::integer>(1) << (2 * m - 1);
        mask_subset = static_cast<dyck::integer>(1) << (set_size - 1);

        MOTZKIN_WORDS.push_back(motzkin_word_t());
        MOTZKIN_WORDS[i].resize(motzkin_word_length);

        for (index_t j = 0; j < motzkin_word_length; j++, mask_subset >>= 1) {
          if (mask_subset & s) {
            MOTZKIN_WORDS[i][j] = j;
          } else {
            if (mask_word & w) {
              stack.push(j);
            } else {
              MOTZKIN_WORDS[i][j]           = stack.top();
              MOTZKIN_WORDS[i][stack.top()] = j;
              stack.pop();
            }
            mask_word >>= 1;
          }
        }
        MOTZKIN_OUTER.push_back(motzkin_word_t());
        MOTZKIN_BOOL.push_back(std::vector<bool>(motzkin_word_length, false));
        for (index_t j = 0; j < set_size; j = MOTZKIN_WORDS[i][j], j++) {
          if (j != MOTZKIN_WORDS[i][j] && MOTZKIN_WORDS[i][j] < set_size) {
            MOTZKIN_OUTER[i].push_back(j);
            MOTZKIN_BOOL[i][j] = true;
          }
        }
        i++;
      }
    }
  }
}

void print_mem_usage(Timer& timer) {
  timer.print();
  std::cout << std::endl;

  double mem = 0;

  for (auto const& x : MOTZKIN_OUTER) {
    mem += (x.size() * sizeof(letter_t));
  }
  size_t len = MOTZKIN_WORDS.at(0).size();

  mem += (MOTZKIN_WORDS.size() * sizeof(letter_t) * len);
  mem += (MOTZKIN_WORDS.size() * sizeof(bool) * len);

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
  std::cout << "Motzkin words use ~ " << mem << suf << std::endl;
  std::cout << "Using " << nr_threads << " / "
            << std::thread::hardware_concurrency() << " threads" << std::endl;
}

void distribute_to_threads_v1(std::vector<std::vector<index_t>>& unprocessed) {
  size_t const nr_motzkin_words = MOTZKIN_WORDS.size();
  size_t const av_load =
      (nr_motzkin_words * (nr_motzkin_words - 1)) / (2 * nr_threads);
  size_t thread_id   = 0;
  size_t thread_load = 0;

  for (size_t i = 0; i < nr_threads; i++) {
    unprocessed.push_back(std::vector<index_t>());
  }

  for (index_t i = 0; i < nr_motzkin_words; i++) {
    unprocessed[thread_id].push_back(i);
    thread_load += nr_motzkin_words - i - 1;
    if (thread_load >= av_load && thread_id != nr_threads - 1) {
      thread_id++;
      thread_load = 0;
    }
  }
}

void distribute_to_threads_v2(std::vector<std::vector<index_t>>& unprocessed) {
  size_t const nr_motzkin_words = MOTZKIN_WORDS.size();
  size_t const av_load =
      (nr_motzkin_words * (nr_motzkin_words - 1)) / (2 * nr_threads);
  std::vector<size_t> thread_load(nr_threads, 0);
  std::vector<size_t> map;
  map.reserve(nr_threads);
  // map from random number to thread_id

  for (size_t i = 0; i < nr_threads; i++) {
    unprocessed.push_back(std::vector<index_t>());
    map.push_back(i);
  }

  std::random_device              rd;
  std::mt19937                    gen(rd());
  std::uniform_int_distribution<> dis(0, nr_threads - 1);
  size_t                          nr_threads_remaining = nr_threads;
  index_t                         i;
  for (i = 0; i < nr_motzkin_words && nr_threads_remaining > 1; i++) {
    size_t thread_id = map[dis(gen)];
    unprocessed[thread_id].push_back(i);
    thread_load[thread_id] += nr_motzkin_words - i - 1;
    if (thread_load[thread_id] >= av_load) {
      nr_threads_remaining--;
      dis = std::uniform_int_distribution<>(0, nr_threads_remaining - 1);
      for (size_t j = thread_id; j < nr_threads - 1; j++) {
        map[j] = map[j + 1];
      }
    }
  }
  size_t thread_id = map[0];
  for (; i < nr_motzkin_words; i++) {
    unprocessed[thread_id].push_back(i);
  }
}

void count_even_rank(size_t                      thread_id,
                     size_t                      nr_motzkin_words,
                     std::vector<index_t> const& unprocessed,
                     size_t&                     nr_idempotents) {
  Timer timer;
  if (verbose) {
    timer.start();
  }
  for (index_t i : unprocessed) {
    nr_idempotents += pow(2, MOTZKIN_OUTER[i].size());
    for (index_t j = i + 1; j < nr_motzkin_words; j++) {
      size_t                         max = 0, cnt = 1;
      motzkin_word_t::const_iterator it = MOTZKIN_OUTER[j].begin();
      do {
        while (*it < max) it++;

        size_t pos  = *it;
        size_t nr_i = 0, nr_j = 1;
        bool   stop = false;

        max = MOTZKIN_WORDS[j][pos];

        if (MOTZKIN_BOOL[i][pos]) {
          nr_i++;
        }

        pos = MOTZKIN_WORDS[j][pos];

        if (pos == MOTZKIN_WORDS[i][pos]) {
          continue;
        }

        pos = MOTZKIN_WORDS[i][pos];

        while (*it != pos) {
          if (MOTZKIN_BOOL[j][pos]) {
            nr_j++;
            if (MOTZKIN_WORDS[j][pos] > max) {
              max = MOTZKIN_WORDS[j][pos];
            }
          } else if (MOTZKIN_BOOL[i][pos]) {
            nr_i++;
          }
          // Check if we reached a fixed point
          if (pos == MOTZKIN_WORDS[j][pos]) {
            stop = true;
            break;
          }
          pos = MOTZKIN_WORDS[j][pos];
          if (pos == MOTZKIN_WORDS[i][pos]) {
            stop = true;
            break;
          }
          pos = MOTZKIN_WORDS[i][pos];
        }
        if (!stop) {
          cnt *= (nr_i * nr_j + 1);
        }
      } while (max < MOTZKIN_OUTER[j].back());
      nr_idempotents += (2 * cnt);
    }
  }
  if (verbose) {
    mtx.lock();
    std::cout << "Thread " << thread_id
              << " is finished, elapsed time = " << timer.string() << std::endl;

    mtx.unlock();
  }
}

// count idempotents of rank 1

void count_odd_rank(size_t                      thread_id,
                    size_t                      nr_motzkin_words,
                    size_t                      deg,
                    std::vector<index_t> const& unprocessed,
                    size_t&                     nr_idempotents) {
  Timer timer;
  if (verbose) {
    timer.start();
  }
  for (index_t const& i : unprocessed) {
    nr_idempotents += pow(2, MOTZKIN_OUTER[i].size());

    for (index_t j = i + 1; j < nr_motzkin_words; j++) {
      // check if there are any idempotents corresponding to the Motzkin words
      // i and j
      size_t pos           = deg;
      bool   no_idempotent = false;
      do {
        if (MOTZKIN_WORDS[i][pos] == pos) {
          no_idempotent = true;
          break;
        }
        pos = MOTZKIN_WORDS[i][pos];
        if (MOTZKIN_WORDS[j][pos] == pos) {
          no_idempotent = true;
          break;
        }
        pos = MOTZKIN_WORDS[j][pos];
      } while (pos != deg);

      if (no_idempotent) {
        continue;
      } else if (MOTZKIN_OUTER[j].empty() || MOTZKIN_OUTER[i].empty()) {
        nr_idempotents += 2;
        continue;
      }

      size_t                         max = 0, cnt = 1;
      motzkin_word_t::const_iterator it = MOTZKIN_OUTER[j].begin();
      do {
        while (*it < max) it++;
        size_t pos  = *it;
        size_t nr_i = (MOTZKIN_BOOL[i][pos] ? 1 : 0);
        size_t nr_j = 1;
        bool   stop = false;

        pos = MOTZKIN_WORDS[j][pos];
        max = pos;

        if (pos != MOTZKIN_WORDS[i][pos]) {
          pos = MOTZKIN_WORDS[i][pos];

          while (*it != pos) {
            if (MOTZKIN_BOOL[j][pos]) {
              nr_j++;
              if (MOTZKIN_WORDS[j][pos] > max) {
                max = MOTZKIN_WORDS[j][pos];
              }
            } else if (MOTZKIN_BOOL[i][pos]) {
              nr_i++;
            }
            // Check if we reached a fixed point
            if (pos == MOTZKIN_WORDS[j][pos]) {
              stop = true;
              break;
            }
            pos = MOTZKIN_WORDS[j][pos];
            if (pos == MOTZKIN_WORDS[i][pos] || pos == deg) {
              stop = true;
              break;
            }
            pos = MOTZKIN_WORDS[i][pos];
          }
          if (!stop) {
            cnt *= (nr_i * nr_j + 1);
          }
        }
      } while (max < MOTZKIN_OUTER[j].back() && it != MOTZKIN_OUTER[j].end());
      nr_idempotents += (2 * cnt);
    }
  }
  if (verbose) {
    mtx.lock();
    std::cout << "Thread " << thread_id
              << " is finished, elapsed time = " << timer.string() << std::endl;

    mtx.unlock();
  }
}

void verify() {
  assert(MOTZKIN_OUTER.size() == MOTZKIN_WORDS.size());
  assert(MOTZKIN_BOOL.size() == MOTZKIN_WORDS.size());

  for (size_t i = 0; i < MOTZKIN_WORDS.size(); i++) {
    assert((size_t) std::count(
               MOTZKIN_BOOL[i].begin(), MOTZKIN_BOOL[i].end(), true)
           == MOTZKIN_OUTER[i].size());
    for (auto j : MOTZKIN_OUTER[i]) {
      assert(MOTZKIN_WORDS[i][j] != j);
    }
  }
}

int main(int argc, char* argv[]) {
  // Not very robust parsing!
  int64_t deg = 0;
  verbose     = false;

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

  index_t n;
  if ((deg / 2) * 2 == deg) {  // deg is even
    n = deg / 2;  // input to dyck, half the length of the returned words
  } else {
    n = (deg - 1) / 2;
  }

  Timer gtimer;
  gtimer.start();

  size_t nr_even_rank = 0;
  size_t nr_odd_rank  = 0;

  {  // Count even rank idempotents
    size_t nr_motzkin_words = nr_motzkin_words_weight_0[deg];

    Timer timer;
    if (verbose) {
      std::cout << "Counting even rank Motzkin idempotents . . ." << std::endl;
      std::cout << "Number of weight 0 Motzkin words is " << nr_motzkin_words
                << std::endl;
      std::cout << "Processing Motzkin words, elapsed time = ";
      timer.start();
    }
    // number of idempotents corresponding to the empty Dyck word
    nr_even_rank += 2 * nr_motzkin_words - 1;
    // don't consider the Motzkin word corresponding to the empty Dyck word
    nr_motzkin_words--;

    if ((deg / 2) * 2 == deg) {
      auto subset_size = [n](size_t m) { return 2 * n - 2 * m; };

      init_motzkin(nr_motzkin_words, 2 * n, 1, n, 2 * n, subset_size);
    } else {
      auto subset_size = [n](size_t m) { return 2 * n - 2 * m + 1; };

      init_motzkin(nr_motzkin_words, 2 * n + 1, 1, n, 2 * n + 1, subset_size);
    }

    if (verbose) {
      print_mem_usage(timer);
    }

    std::vector<std::vector<index_t>> unprocessed;
    distribute_to_threads_v1(unprocessed);

    std::vector<size_t>      nr_idempotents(nr_threads, 0);
    std::vector<std::thread> threads;

    for (size_t i = 0; i < nr_threads; i++) {
      threads.push_back(std::thread(count_even_rank,
                                    i,
                                    nr_motzkin_words,
                                    std::ref(unprocessed[i]),
                                    std::ref(nr_idempotents[i])));
    }

    // corresponds to empty Dyck words and whole set as subset
    for (size_t i = 0; i < nr_threads; i++) {
      threads[i].join();
      nr_even_rank += nr_idempotents[i];
    }

    if (verbose) {
      std::cout << "There are " << nr_even_rank << " even rank idempotents, ";
      std::cout << "elapsed time = ";
      timer.print();
      std::cout << std::endl;
    }
  }
  {  // Count odd rank idempotents
    size_t nr_motzkin_words = nr_motzkin_words_weight_1[deg];
    Timer  timer;

    if (verbose) {
      std::cout << "Counting odd rank Motzkin idempotents . . ." << std::endl;
      std::cout << "Number of weight 1 Motzkin words is " << nr_motzkin_words
                << std::endl;
      std::cout << "Processing Motzkin words, elapsed time = ";
      timer.start();
    }
    if ((deg / 2) * 2 == deg) {
      auto subset_size = [n](size_t m) { return 2 * n - 2 * m + 1; };

      init_motzkin(nr_motzkin_words, 2 * n + 1, 1, n, 2 * n, subset_size);
    } else {
      auto subset_size = [n](size_t m) { return 2 * n - 2 * m + 2; };

      init_motzkin(
          nr_motzkin_words, 2 * n + 2, 1, n + 1, 2 * n + 1, subset_size);
    }

    if (verbose) {
      print_mem_usage(timer);
    }

    std::vector<std::vector<index_t>> unprocessed;
    distribute_to_threads_v1(unprocessed);

    std::vector<size_t>      nr_idempotents(nr_threads, 0);
    std::vector<std::thread> threads;

    for (size_t i = 0; i < nr_threads; i++) {
      threads.push_back(std::thread(count_odd_rank,
                                    i,
                                    nr_motzkin_words,
                                    deg,
                                    std::ref(unprocessed[i]),
                                    std::ref(nr_idempotents[i])));
    }

    // corresponds to empty Dyck words and whole set as subset
    for (size_t i = 0; i < nr_threads; i++) {
      threads[i].join();
      nr_odd_rank += nr_idempotents[i];
    }
    if (verbose) {
      std::cout << "There are " << nr_odd_rank << " odd rank idempotents, ";
      std::cout << "elapsed time = ";
      timer.print();
      std::cout << std::endl;
    }
  }

  if (verbose) {
    std::cout << "Total elapsed time = ";
    gtimer.print();
    std::cout << std::endl;
  }

  std::cout << nr_even_rank + nr_odd_rank << std::endl;
  exit(0);
}
