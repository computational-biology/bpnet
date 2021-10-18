NAME	
    BPNet


SYNOPSIS
		bpnet.linux [OPTIONS] [FILE...]


DESCRIPTION
	This program is for calculating the overlap between 
	two bases in a nucluic acid for RNA. The output is
	stored in a .rob file with the same name of the cor
	file name. 


INSTLLATION
      	EASY-METHOD
		Download the binary executable 'bpnet.linux' and place it to some suitable folder.
		Done.

	HARD-METHOD
		If that does not work for you for any reason, then download the src folder 
		and compile the same following the instructions given below.

		COMPILE:
			cd src
			make clean
			make
		
		Then move the binary executable 'bpnet.linux' to your suitable path.

SETUP
	BPNet uses several files for its work. These are kept in 'sysfiles' directory.
	copy all of them and paste them in a folder. We prefer the following.

		/usr/local/bin

PATH-SETTING
	BPNet reads system files through an environment variable called NUCLEIC_ACID_DIR. 
	Suppose you have put all the files of 'sysfiles' directory in the /usr/local/bin, 
	then put the following line in your '.bashrc' file.

		export NUCLEIC_ACID_DIR=/usr/local/bin/
	
PARALLEL-MODE
	BPNet can run in parallel mode. But for this OpenMP library should be installed
	in your machine.
		PARALLEL-MODE COMPILE
			cd src
			make clean
			make -f Makefileprll
	
		PARALLEL-MODE SETUP
			In parallel mode, you have to do all the setup and path settings
			of the serial mode. In addition you have to set the OpenMP
			environment variable. If our machine has 4 cores and each having
			4 hyperthreads, i.e. 2x4 = 8 cpu. Then to achieve maximum from
			the system, set the path as,

				export OMP_NUM_THREADS=8


RUN
	To run the program, go to the directory where the structure files are stored.
	Then run as follows.
		
		bpnet.linux xxxx.cif

INPUT
	BPNet accepts mmCIF and PDB file formats.

OUTPUT
	A full list of output files are listed below.

		xxxx_rna.pdb 
			This file stores the RNA part only in PDB format.
		xxxx.out
			This file stores all the base-pair information.
		xxxx.adj
			It stores the adjacency information.
		xxxx.fasta
			This file stores the nucleic acid primary sequence in fasta format.
		xxxx.dbn
			This file is the base pair information in dot-bracket notation.
		xxxx.bpseq
			This file is the base pair information in bpseq notation.
		xxxx.dat
			This file stores the secondary structure information in fasta format. 
			The sequence here will look like, HHHHCCLLL, in this way. Here H means 
			the corresponding base is in helix, C represent that the residue is in 
			Coil, L Stands for Loop etc.
		xxxx.hlx 
			This files stores the helix, loop and pseudo helix information.
		xxxx_rna.pml 
			This one is to visualize the networks in RNA part only.
		xxxx_cif.pml
			This file is to visualize the same on the entire mmCIF file's context.
			This program is an applet file ready for visualization using VARNA [Ref??] software.
			This file stores all helix information.
			A postscript file is generated to visualize the base base interaction.
			The program also generates two secondary structure file namely dbn and bpseq file.
		xxxx.ps
			A postscript file to show the contact-map distribution.
		xxxx_helix.pml
			To visualize the networks in the context of other functional units like base pairs,
			helix, loops etc.
		xxxx.html
			An applet file is generated to visualize the secondary structures in a better way.
			The applet is compatible with VARNA (http://varna.lri.fr). 

RUNTIME OPTIONS
	BPNet uses several command line options to tune the result. The lollowing are the full lists 
	of the same.
	
	GENERAL OPTIONS

		--help
			Shows a small help on terminal.

		--genhelp
			Generates this help.md file in the current directory.

		--version
			Prints the version of the program.

		-nettype=[basepair/contact]
			This option indicates base pair or contact based network.
			Default is 'contact' based network.

			EXAMPLE	
				bpnet.linux  xxxx.cif -nettype=basepair

		-netzise=from-to    or   -netsize=size
			This option generates the network of the mentioned size.
			default is 3-999999.

			EXAMPLE
				bpnet.linux  xxxx.cif -netsize=4-5
				bpnet.linux  xxxx.cif -netsize=3

		-exdeg=number
			This option generates all the network which contains at least
			one vertex with degree 'count'. 'exdeg' means 'exists degree'.

			EXAMPLE 
				bpnet.linux  xxxx.cif -exdeg=3

		-numexdeg=number
			This one considers those networks where the total nodes with specific 
			degree is required.

			EXAMPLE
				bpnet.linux xxxx.cif -exdeg=3 -numexdeg=2

		-cycles=number
			This option generates the networks where there is at least two cycles.

			EXAMPLE
				bpnet.linux xxxx.cif -cycles=2

		-wttype=[c1p-c1p]
			Default BPNet considers E-val (Ref: https://doi.org/10.1080/07391102.2006.10507108)
			for base pairs and surface overlap for contact based network. But these can be changed
			to C1'-C1' distance. 
	
			EXAMPLE
				bpnet.linux xxxx.cif -wttype=c1p-c1p
				bpnet.linux xxxx.cif -wttype=c1p-c1p -nettype= basepair

		-wtcutoff=value
			Default the contact based network considers any bases whose surface
			contact greater than 0. But the user bay change it to a desired positive
			value. In that case all contacts greater ot equal to that value will be
			considered for contact based network.

			EXAMPLE
				bpnet.linux xxxx.cif -wtcutoff=10.0


	BASEPAIR RELATED OPTIONS
		BPNet accepts several options for base pair generations.
		
		-hbdist=value
			BPNet program uses 3.8A as the default distance of donor-acceptor
			atoms for Hydrogen bond. User may change the same using this option.

			EXAMPLE
				bpnet.linux xxxx.cif -bhdist=4.2
		
		-sugmed=[true/false]
			BPNet by default considers O2' atom of sugar-phosphate backbone for sugar
			mediated base pairs. If the user does not want that, then they can set it
			to 'false'

			EXAMPLE
				bpnet.linux xxxx.cif -sugmed=false

		-chmed=[true/false]
			By default BPNet considers C-H...O/N mediated H-bonds for base pairs. It one
			wants to exclude them then they can choose it as false.

			EXAMPLE
				bpnet.linux xxxx.cif -chgmed=false

		-hetatm=[true/false]
			By default BPNet considers all modified nucleobases. If the user wants
			to exclude them, then they can do so by this option.
				
			EXAMPLE
				bpnet.linux xxxx.cif -hetatm=false

		-chain=name
			If one supply the chain name, the program computes network on the
			mentioned chain.

			CAUTION !!!
				If a chain is selected then the program only sees that chain.
				So, if a base of the given chain has a network with the bases 
				of the other chain, the program will not report that.

			EXAMPLE
				bpnet.linux xxxx.cif -chain=A
		

ABBREVIATIONS
	
	BASE PAIR EDGE

		W - Watson-Crick edge (Capital W).
		H - Hoogsteen edge (Capital H).
		S - Sugar edge (Capital S).
    		w - Watson-Crick edge with one or more C-H...O/N type of hydrogen bond (Small w).
    		h - Hoogsteen edge with one or more C-H...O/N type of hydrogen bond (Small h).
    		s - Sugar edge with one or more C-H...O/N type of hydrogen bond (Small s).
    		+ - Protonated Watson-Crick edge.
    		z - Protonated Sugar edge.
    		g - Protonated Hoogsteen edge (rarely found though).
	
	ORIENTATION
    		C - Cis Orientation.
    		T - Trans Orientation.

	PAIR TYPE

    		BP - Normal base pair.
    		TP - Tartiary pair.
    		BF - Bifurcated pair (Follow our paper BPFIND(2006) )

BUG REPORT
		Email bug reports to the bug-reporting address 
		⟨roy.parthajit@gmail.com⟩  or 
		⟨dhananjay.bhattacharyya@saha.ac.in⟩

				





