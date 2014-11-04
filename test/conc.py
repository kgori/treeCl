import treeCl

c = treeCl.Collection(input_dir="/Users/kgori/scratch/simtest4", file_format="phylip")
print c[5].get_distances()
print c[5].chkdst()
print c[5].get_bionj_tree()
#conc = treeCl.Concatenation(c,[5,6,7,8,9])
#part = conc.qfile(protein_model='LGX')
#p = conc.alignment.pll_get_instance(part, c[5].tree, 6)
