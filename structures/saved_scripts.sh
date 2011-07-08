

for f in *.pdb; do cat $f | str cent | str rup $num > ${f%%.pdb}_$num.pdb; done; for f in *$num*pdb; do molauto  $f | mol lab_r $num $((num+5)) | molscript > ${f%%.pdb}.ps; done; open *.ps
