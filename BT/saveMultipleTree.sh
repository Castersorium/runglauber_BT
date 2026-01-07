# for alpha in 2.0 3.0 4.0; do
#   root -l -b -q "saveTreeAA.C(\"Au197pnHFB14\",\"Au197pnHFB14\",\"200\",${alpha},0,5)"
# done

for alpha in  2.0 3.0 3.25 3.5 3.75 4.0 4.25 4.5 ; do
  root -l -b -q "saveTreeAA.C(\"Pbpnrw\",\"Pbpnrw\",\"17p3\",${alpha},0,5)"
  echo "Complete Saving alpha = ${alpha}"
done

for alpha in  2.0 3.0 3.25 3.5 3.75 4.0 4.25 4.5 ; do
  root -l -b -q "readTree.C(\"Pbpnrw\",\"Pbpnrw\",\"17p3\",${alpha},0,5)"
  echo "Complete Reading alpha = ${alpha}"
done