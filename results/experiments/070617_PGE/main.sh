#mysql -u anonymous -h ensembldb.ensembl.org -e "show databases;" > databases.txt
#mysql -sN -u anonymous -h ensembldb.ensembl.org -e   "use homo_sapiens_core_75_37; SELECT gene.stable_id, seq_region.name, gene.seq_region_start, gene.seq_region_end, gene.description FROM gene, seq_region WHERE gene.seq_region_id=seq_region.seq_region_id AND seq_region.name IN ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT')" > ensembl75.txt


perl pge.pl -q parental_vs_resistant.txt -r ensembl75 > ./parental_vs_resistant/parental_vs_resistant_results.txt
perl pge.pl -q parental_vs_resistant_up.txt -r ensembl75 > ./parental_vs_resistant/parental_vs_resistant_results_up.txt
perl pge.pl -q parental_vs_resistant_down.txt -r ensembl75 > ./parental_vs_resistant/parental_vs_resistant_results_down.txt


perl pge.pl -q parental_vs_resistant_Met+Hcy-.txt -r ensembl75 > ./parental_vs_resistant_Met+Hcy-/parental_vs_resistant_Met+Hcy-_results.txt
perl pge.pl -q parental_vs_resistant_Met+Hcy-_up.txt -r ensembl75 > ./parental_vs_resistant_Met+Hcy-/parental_vs_resistant_Met+Hcy-_results_up.txt
perl pge.pl -q parental_vs_resistant_Met+Hcy-_down.txt -r ensembl75 > ./parental_vs_resistant_Met+Hcy-/parental_vs_resistant_Met+Hcy-_results_down.txt


perl pge.pl -q parental_vs_resistant_Met-Hcy+.txt -r ensembl75 > ./parental_vs_resistant_Met-Hcy+/parental_vs_resistant_Met-Hcy+_results.txt
perl pge.pl -q parental_vs_resistant_Met-Hcy+_up.txt -r ensembl75 > ./parental_vs_resistant_Met-Hcy+/parental_vs_resistant_Met-Hcy+_results_up.txt
perl pge.pl -q parental_vs_resistant_Met-Hcy+_down.txt -r ensembl75 > ./parental_vs_resistant_Met-Hcy+/parental_vs_resistant_Met-Hcy+_results_down.txt
