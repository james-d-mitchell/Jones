The number of idempotents in a Jones monoid
=====

If `n` is any positive integer, then
```jones n```
returns the number of idempotents in the Jones monoid of degree `n`. If you
want more information about what is going on, then do 
```jones -v n```

`jones` is a multi-threaded C++ program.  You can alter then number of threads
by changing the variable `nr_threads` in line 44 of `jones.cc` to any positive
integer.  `jones` will then use `nr_threads + 1`  threads, regardless of the
maximum number of threads which your hardware supports.

### Compiling

`jones` requires the code from the git repository:

   [https://github.com/cassioneri/Dyck](https://github.com/cassioneri/Dyck)

to be put inside a subdirectory (of the directory where `jones` will be compiled)
named `Dyck`.

Compile with:

```g++ -O3 -std=c++11 -Wall -Wextra -pedantic -o jones jones.cc```

On linux you might have to do:

```g++ -O3 -pthread -std=c++11 -Wall -Wextra -pedantic -o jones jones.cc```

### Testing

You can test if `jones` is working properly by running the `test.sh` script. If nothing is returned/displayed, then everything is working.

Enjoy!

Copyright (C) 2016 James D. Mitchell

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.

