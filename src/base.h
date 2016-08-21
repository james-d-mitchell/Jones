/*******************************************************************************

 Copyright (C) 2016 James D. Mitchell

 This work is licensed under a Creative Commons Attribution-ShareAlike 4.0
 International License. See
 http://creativecommons.org/licenses/by-sa/4.0/

*******************************************************************************/

#ifndef BASE_H_
#define BASE_H_

#include <bitset>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "timer.h"
#include "Dyck/dyck.h"

// Typedefs

typedef uint_fast8_t          letter_t;
typedef std::vector<letter_t> dyck_vec_t;
typedef size_t                dyck_index_t;

static size_t const catalan_numbers[] = {1, 1, 2, 5, 14, 42, 132, 429, 1430,
  4862, 16796, 58786, 208012, 742900, 2674440, 9694845, 35357670, 129644790,
  477638700, 1767263190, 6564120420, 24466267020, 91482563640, 343059613650,
  1289904147324, 4861946401452, 18367353072152, 69533550916004,
  263747951750360, 1002242216651368, 3814986502092304};

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

void print_binary(dyck::integer val) {
  std::cout << std::bitset<sizeof(dyck::integer) * 8>(val).to_string()
            << std::endl;
}

std::string string_mem(double mem) {
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
  return std::to_string(mem) + suf;
}

void parse_args(int argc, char* argv[], bool& verbose, size_t& deg) {
  // Not very robust parsing!
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
          argv[0] << ": must be an integer in [1, 40]" << std::endl;
        exit(-1);
      }
    }
  }
}


// reverse bits in w
dyck::integer reverse(dyck::integer w, size_t dyck_word_length) {
  size_t        nr_bits = sizeof(dyck::integer) * 8;
  dyck::integer a       = 0;
  unsigned char b;

  while (nr_bits > 0) {
    // reverse the next byte of w
    b = w;
    b = ~(((b * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32);
    nr_bits -= 8;
    a |= ((dyck::integer)b << nr_bits);
    w >>= 8;
  }
  return (a >> (sizeof(dyck::integer) * 8 - dyck_word_length));
}

#endif  // BASE_H_
