#!/bin/bash

energy="17p3_2011"
projectile="Pbpnrw"
target="Pbpnrw"

# alpha 扫描
alphas=(3.0 4.0 5.0)

# 中心度区间（百分比）
centMins=(0    5     12.5  23.5  33.5)
centMaxs=(5    12.5  23.5  33.5  43.5)

for alpha in "${alphas[@]}"; do
  for i in ${!centMins[@]}; do
    cmin=${centMins[$i]}
    cmax=${centMaxs[$i]}

    root -l -b -q \
      "saveTreeAA.C(\"${projectile}\",\"${target}\",\"${energy}\",${alpha},${cmin},${cmax})"
    
    root -l -b -q \
      "readTree.C(\"${projectile}\",\"${target}\",\"${energy}\",${alpha},${cmin},${cmax})"

    echo "Complete: alpha=${alpha}, cent=${cmin}-${cmax}%"
  done
done
