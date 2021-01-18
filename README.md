This program is for finding the base-pair networks as well as overlap networks In DNA/RNAs.
 It assumes a cif, pdb or out file as its input.

 The switches are:
 
     -netsize=5 (prints only network of size 5), Default is 3
     
     -netsize=5-10 (prints only network of size 5 to 10), Default is 3-30
     
     -exdeg=3  (Prints ann networks with at least a vertex with degree 3
                (Default is 2)
                
     -cycles=2 (Computes the number of cycles, Main cycles), Default 0
     
     -ovlpcutoff=20.0 for non-base pair overlap cutoff value
     
     -adj=true for creation of adjacency matrix. Default false


 It generates a .edge file for edge lists, .adj file for adjacency
 matrix for every file. It also generates a pairchain.net for all files.
 It works for multiple files, like, 
 
     bpnetl.linux *.cif -netsize=4 -exdeg=3
 
 Installation: Follow the installation_guide.txt
 
 Compilation:  Follow the installation_guide.txt
 
 run:    bpnet.linux [optional switches] 1ehz.cif 


