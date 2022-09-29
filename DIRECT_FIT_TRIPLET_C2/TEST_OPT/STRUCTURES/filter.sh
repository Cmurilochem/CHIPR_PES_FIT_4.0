#!/bin/bash

nS=$(grep -n -i "Molden Format" ../min.molden | tail -n 1 | gawk '{print $4}')

nSm=$((nS-1))

for i in `seq -w 1 ${nSm}`; do
    gawk -v nst=${i} '{if ( /Molden Format/ && $3==nst ){print $0;getline;while($0 !~ /Molden Format/){print $0; getline}}}' ../min.molden > min_molden_$i.res
done

last=$(grep -n -i "Molden Format" ../min.molden | cut -d: -f1 | tail -n 1)
total=$(wc -l ../min.molden | cut -d\  -f1)

tail -n $((total-last+1)) ../min.molden > min_molden_${nS}.res

exit
