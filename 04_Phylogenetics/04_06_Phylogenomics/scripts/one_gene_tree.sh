grep $gene$ data/part_1/mitochondrions_renamed_final/*mt_prots.fa --no-filename | cut -c2- > data/part_1/gene_names/$(basename "$gene")_names.txt
for file in data/part_1/mitochondrions_renamed_final/*.mt_prots.fa; do bbmap/filterbyname.sh in=$file out=data/part_1/mt_genes/sep/$(basename "$file" .mt_prots.fa).$gene.aa.fa include=t names=data/part_1/gene_names/$(basename "$gene")_names.txt overwrite=true ignorejunk=true; done
cat data/part_1/mt_genes/sep/*$gene.aa.fa > data/part_1/mt_genes/merged/$gene.aa.fa
mafft --auto data/part_1/mt_genes/merged/$gene.aa.fa > data/part_1/mt_aligns/$gene.aa.aln
