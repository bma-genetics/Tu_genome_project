
1. List collinear genes identified between Tu3, Ta3B, Brachypodium, rice and sorguhum within duplication blocks containing > 5 genes between Tu3 and Ta3B

perl format_collinear.pl  ====> Tu_ref_collineariy.out/Ta_ref_collinearity.out


2. Identify gene insertions and deletions on the basis of Tu/Ta_ref_collineariy

id_InDel.pl ====> Ta_ref_collinearity_indel.out/Tu_ref_collinearity_indel.out


3. Count insertions and deletions for Tu and Ta genes

count_InDel.pl ===> Ta_insertion 648; Tu_insertion 354; Tu_deletion 393; Ta_deletion 213
