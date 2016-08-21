Counting idempotents in a monoid of partitions
=====

The code in this project can be used to count idempotents in the Jones, Motzkin, and Kauffman monoids, and is based on [Idempotent statistics of the Motzkin and Jones monoids](http://arxiv.org/abs/1507.04838) by Igor Dolinka, James East, Athanasios Evangelou, Desmond FitzGerald, Nicholas Ham, James Hyde, and Nicholas Loughlin. 

###Installation
In your terminal of choice, perform the following steps to install this project.

0. Clone this repository by doing 

        git clone https://github.com/james-d-mitchell/Jones
      
    This will create a directory called `Jones` in whatever directory you were in to begin with.
1. Type `cd Jones` to change directory into the `Jones` folder created in Step 1.
2. Type 

        git clone https://github.com/cassioneri/Dyck
   
   This will create a directory called `Dyck` inside the directory `Jones`. 
3. Type `make` which will create executables for the three programs in this project `jones`, `motzkin` and `kauffman`. This project is written in C++11, and so you will require a C++ compiler compatible witht this standard.

You might want to test that everything worked by doing `make test` which should display 

    tst/jones.sh
    tst/motzkin.sh
    tst/kauffman.sh
    
but nothing further.

### Counting idempotents

If `n` is any positive integer, then
`jones n`
returns the number of idempotents in the Jones monoid of degree `n`. If you
want more information about what is going on, then do 
`jones -v n`
`motzkin` and `kaufmann` can be used in the same way.

`jones`, `motzkin`, and `kauffman` are a multi-threaded C++ programs. By default the number of threads used is one less than the maximum supported by the hardware. 

You can alter the number of threads by changing the variable `nr_threads` in the files 
`src/jones.cc`, `src/kauffman.cc`, and/or `src/motkzin.cc`, and rerunning `make`.
Each of the programs will then use `nr_threads + 1` threads, regardless of the
maximum number of threads which your hardware supports.

Enjoy!

Copyright (C) 2016 James D. Mitchell

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.

