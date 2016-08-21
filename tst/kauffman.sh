#!/bin/bash
set -e
if [ ! -f ./kauffman ]; then
  echo "kauffman executable not found, please build it!"
  exit 1
fi

if [ -f tst/results ]; then
  rm -f tst/results
fi

for i in {1..16}
do
  ./kauffman $i >> tst/results
done

diff tst/results tst/expected-kauffman

if [ -f tst/results ]; then
  rm -f tst/results
fi
