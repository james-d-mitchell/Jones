#!/bin/bash
set -e
if [ ! -f ./jones ]; then
  echo "jones executable not found, please build it!"
  exit 1
fi

if [ -f tst/results ]; then
  rm -f tst/results
fi

for i in {1..18}
do
  ./jones $i >> tst/results
done

diff tst/results tst/expected-jones

if [ -f tst/results ]; then
  rm -f tst/results
fi

