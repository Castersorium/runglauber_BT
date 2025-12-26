for alpha in 1.0 2.0 3.0 4.0 5.0; do
  root -l -b -q "saveTreeAA.C(\"Au197pnHFB14\",\"Au197pnHFB14\",\"62p4\",${alpha},0,5)"
done

for alpha in 1.0 2.0 3.0 4.0 5.0; do
  root -l -b -q "saveTreeAA.C(\"Au197pnHFB14\",\"Au197pnHFB14\",\"62p4\",${alpha},70,80)"
done
