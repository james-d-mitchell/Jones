#!/bin/bash
if [ -f results ]; then
  rm -f results
fi

for i in {1..18}
do
  ./jones $i >> results
done
diff results expected
