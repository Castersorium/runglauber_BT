#!/bin/bash

# You may:
# sed -i 's/\r$//' saveMultipleTree.sh
# first before run this script

set -e

# ================= User control =================
ENERGY_MODE=( "62p4_STAR" )   # 17p3_NA49 | 62p4_BRAH | 200_BRAH
ALPHAS=(3.0 4.0 5.0)        # just edit here
# ===============================================

case "$ENERGY_MODE" in
  17p3_NA49)
    energy="17p3_1999"
    projectile="Pbpnrw"
    target="Pbpnrw"
    centMins=(0)
    centMaxs=(5)
    ;;
  62p4_BRAH)
    energy="62p4_BRAH"
    projectile="Au197pnHFB14"
    target="Au197pnHFB14"
    centMins=(0)
    centMaxs=(10)
    ;;
  62p4_STAR)
    energy="62p4_STAR"
    projectile="Au197pnHFB14"
    target="Au197pnHFB14"
    centMins=(0)
    centMaxs=(80)
    ;;
  200_BRAH)
    energy="200_BRAH"
    projectile="Au197pnHFB14"
    target="Au197pnHFB14"
    centMins=(0)
    centMaxs=(5)
    ;;
  200_STAR)
    energy="200_STAR"
    projectile="Au197pnHFB14"
    target="Au197pnHFB14"
    centMins=(0)
    centMaxs=(80)
    ;;
  *)
    echo "Unknown ENERGY_MODE: $ENERGY_MODE"
    exit 1
    ;;
esac


for alpha in "${ALPHAS[@]}"; do
  for i in "${!centMins[@]}"; do
    cmin=${centMins[$i]}
    cmax=${centMaxs[$i]}

    echo "Running: ${energy}, alpha=${alpha}, cent=${cmin}-${cmax}%"

    root -l -b -q \
      "saveTreeAA.C(\"${projectile}\",\"${target}\",\"${energy}\",${alpha},${cmin},${cmax})"

    root -l -b -q \
      "readTree.C(\"${projectile}\",\"${target}\",\"${energy}\",${alpha},${cmin},${cmax})"

    echo "Done."
    echo "-------------------------------------------"
  done
done
