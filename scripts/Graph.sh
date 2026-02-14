# run the Haddock
for i in $(ls data/); do python scripts/run_pipeline.py --complexes $i; done

# get the results
for i in $( ls results); do ID=$(echo $i| sed 's/run_//;s/_Epi//'); for ii in $(ls results/$i/2_flexref/*.gz); do NAME=$(echo $ii | awk -F'/' '{print $4}' |sed "s/flexref/$ID/");cp $ii PDBs/$NAME; done ;done
gunzip PDBs/*
for i in $( ls results); do ID=$(echo $i| sed 's/run_//;s/_Epi//'); for ii in $(ls results/$i/2_flexref/*.pdb); do NAME=$(echo $ii | awk -F'/' '{print $4}' |sed "s/flexref/$ID/");cp $ii PDBs/$NAME; done ;done

# get the reference
cp  results/*/complex_merged_*.pdb  PDBs

cp PDBs/complex_merged_* PDB_Epi

# convert into graph
python scripts/flexref_to_graphs.py --runs PDBs
