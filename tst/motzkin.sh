#!/bin/bash
set -e
if [ ! -f ./motzkin ]; then
  echo "motzkin executable not found, please build it!"
  exit 1
fi

if [ -f tst/results ]; then
  rm -f tst/results
fi

for i in {1..11}
do
  ./motzkin $i >> tst/results
done

diff tst/results tst/expected-motzkin

if [ -f tst/results ]; then
  rm -f tst/results
fi
