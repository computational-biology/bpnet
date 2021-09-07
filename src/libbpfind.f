C
C Released version in June 2021, version 2.1.17
C
C
! Copyright 2019 Saha Institute of Nuclear Physics, Kolkata, INDIA
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
! http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the License is distributed on
! an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
! KIND, either express or implied. See the License for the
! specific language governing permissions and limitations
! under the License.

C  Finding base pairing information from a PDB file
C
C  The method of this program was described in the following publication:
C  
C  J. Das, S. Mukherjee, A. Mitra and D. Bhattacharyya (2006) Non-Canonical Base 
C  Pairs and Higher Order Structures in Nucleic Acids: Crystal Structure 
C  Database Analysis J. Biomol. Struct. Dynam. 24, 149-161.
C
C
C  We are thankful to Arvind Marathe, MBU, IISc, Bangalore and Purshotam
C  Sharma, IIIT, Hyderabad for valuable bug-reports. Additional members 
C  contributing to the development are Sukanya Halder and Parijat Majumdar,
C  Debasish Mukherjee, SINP, Kolkata, Parthajit Roy of Burdwan University, India
C  
C  Last updated Jan. 2019
C  Last updated on Jan, 2021 by Parthajit Roy on this lib module
C


       subroutine callbpfindc(cif , accn, ht, hd, 
     1 hdval, ang, angval, ch, sg, 
     2 cor, eval, chain, chainval, 
     3 nmr, nmrval)bind(c, name ='callbpfindc')
       use iso_c_binding, only: c_char , c_null_char
C       character(kind=c_char), dimension(*), intent(in)    ::  argv
C       character(kind=c_char), intent(in) ::  a1(20), a2(20)
       character(kind=c_char,len=1),dimension(512),intent(in) ::cif
       character(kind=c_char,len=1),dimension(512),intent(in) ::accn
       character(kind=c_char,len=1),dimension(512),intent(in) ::ht
       character(kind=c_char,len=1),dimension(512),intent(in) ::hd
       character(kind=c_char,len=1),dimension(512),intent(in) ::hdval
       character(kind=c_char,len=1),dimension(512),intent(in) ::ang
       character(kind=c_char,len=1),dimension(512),intent(in) ::angval
       character(kind=c_char,len=1),dimension(512),intent(in) ::ch
       character(kind=c_char,len=1),dimension(512),intent(in) ::sg
       character(kind=c_char,len=1),dimension(512),intent(in) ::cor
       character(kind=c_char,len=1),dimension(512),intent(in) ::eval
       character(kind=c_char,len=1),dimension(512),intent(in) ::chain
       character(kind=c_char,len=1),dimension(512),intent(in) ::chainval
       character(kind=c_char,len=1),dimension(512),intent(in) ::nmr
       character(kind=c_char,len=1),dimension(512),intent(in) ::nmrval
        ! here I have a string with fixed length
        character*512 val1(50)
        character (len=512) ::params

        integer :: i, narg
         
        narg = 0


C        params = "";  
C        loop_argc: do i=1, 20
C        if ( argc(i) == c_null_char ) then
C           exit loop_argc
C        else
C           params (i:i) = argc(i)
C        end if
C        end do loop_argc
C        write (*,*) "narg=",params
C        read(params,*) narg


        params = "";  
        loop_cif: do i=1, 512
        if ( cif(i) == c_null_char ) then
              exit loop_cif
        else
              params (i:i) = cif(i)
        end if
        end do loop_cif


        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if

        params = "";  
        loop_accn: do i=1, 512
        if ( accn(i) == c_null_char ) then
              exit loop_accn
        else
              params (i:i) = accn(i)
        end if
        end do loop_accn
         write(*,*) "       BASE-PAIRS COMPUTATION STARTS"

c        write(*,*) params

        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if



        params = "";  
        loop_ht: do i=1, 512
        if ( ht(i) == c_null_char ) then
              exit loop_ht
        else
              params (i:i) = ht(i)
        end if
        end do loop_ht

c        write(*,*) params

        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if





        params = "";  
        loop_hd: do i=1, 512
        if ( hd(i) == c_null_char ) then
              exit loop_hd
        else
              params (i:i) = hd(i)
        end if
        end do loop_hd

c        write(*,*) params
        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if


        params = "";  
        loop_hdval: do i=1, 512
        if ( hdval(i) == c_null_char ) then
              exit loop_hdval
        else
              params (i:i) = hdval(i)
        end if
        end do loop_hdval

c        write(*,*) params
        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if

        params = "";  
        loop_ang: do i=1, 512
        if ( ang(i) == c_null_char ) then
              exit loop_ang
        else
              params (i:i) = ang(i)
        end if
        end do loop_ang

c        write(*,*) params

        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if

        params = "";  
        loop_angval: do i=1, 512
        if ( angval(i) == c_null_char ) then
              exit loop_angval
        else
              params (i:i) = angval(i)
        end if
        end do loop_angval

c       write(*,*) params

        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if




        params = "";  
        loop_ch: do i=1, 512
        if ( ch(i) == c_null_char ) then
              exit loop_ch
        else
              params (i:i) = ch(i)
        end if
        end do loop_ch

c        write(*,*) params
        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if





        params = "";  
        loop_sg: do i=1, 512
        if ( sg(i) == c_null_char ) then
              exit loop_sg
        else
              params (i:i) = sg(i)
        end if
        end do loop_sg

c       write(*,*) params
        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if

C       Other options are

c        narg = narg + 1
c        val1(narg) = "-NODAT"


        params = "";  
        loop_cor: do i=1, 512
        if ( cor(i) == c_null_char ) then
              exit loop_cor
        else
              params (i:i) = cor(i)
        end if
        end do loop_cor

c       write(*,*) params
        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if





        params = "";  
        loop_eval: do i=1, 512
        if ( eval(i) == c_null_char ) then
              exit loop_eval
        else
              params (i:i) = eval(i)
        end if
        end do loop_eval


        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if

        params = "";  
        loop_chain: do i=1, 512
        if ( chain(i) == c_null_char ) then
              exit loop_chain
        else
              params (i:i) = chain(i)
        end if
        end do loop_chain

c        write(*,*) params

        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if

        params = "";  
        loop_chainval: do i=1, 512
        if ( chainval(i) == c_null_char ) then
              exit loop_chainval
        else
              params (i:i) = chainval(i)
        end if
        end do loop_chainval

c       write(*,*) params

        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if



        params = "";  
        loop_nmr: do i=1, 512
        if ( nmr(i) == c_null_char ) then
              exit loop_nmr
        else
              params (i:i) = nmr(i)
        end if
        end do loop_nmr

c        write(*,*) params

        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if

        params = "";  
        loop_nmrval: do i=1, 512
        if ( nmrval(i) == c_null_char ) then
              exit loop_nmrval
        else
              params (i:i) = nmrval(i)
        end if
        end do loop_nmrval

c       write(*,*) params

        if(params.ne.'-dummyval') then
              narg = narg + 1
            val1(narg) = params
        end if











        narg = narg + 1
        val1(narg) = "-NONUP"



        narg = narg + 1
        val1(narg) = "-NOHLX"



        narg = narg + 1
        val1(narg) = "-NOCSV"



        narg = narg + 1
        val1(narg) = "-FASTA"




c        narg = narg + 1
c        val1(narg) = "-c1"

        







        call bpfind(narg, val1)

        end subroutine callbpfindc

       subroutine bpfind(narg,val1)
       EXTERNAL VERIFY
       DIMENSION NUCL(96),ires(100000),
     1  hbdist(100000),hbangl(100000),localinfo(30),
     2  bpinfo(30),tpinfo(10000),val2pr(5),
     3  energylc(30),anglinfo(10000),bptype(30)
       
       integer atmN,base,ATGC,KS,basef,NN2,localinfo,nopair,prd,
     1         NN1,NN3,NN4,presd,bsinfo,localinfb,nores,locp(30),
     2         chohb,hetatm,sugarbp,nforce,oligo,nnfpr(40)
       real cutoff,energy,energylc,calenrg,varenergy,cutang,occ
       real anglinfo,varangl,yangle,anglinfb,bf,xa,ya,za,ocp
	logical raretauto
       character*5 dummynm
       CHARACTER*4 NAMED,atmnam,ATOM,chaind,pcd,chainid,pchaind,chan2,
     1   molid
       CHARACTER*3 NUCL,resn,resd,locf(30),adevar,guavar,
     1   cytvar,uravar
       Character*3 feature,bpinfo,bs,bptype,bpinfb,bptypeb,resid,ans(2)
       character*1 type,tpinfo,tp,tpinfb,pos(100000),tpod(30),tmppos,
     1 posins,chan,loct(30),strseq,prntins,inscode,allins,
     2 tmpins
       character*2 bptp(30),locb(30)
       Character*80 filenm,nmpass,title,LINE,ttl(50),aline,filen2,
     1 storedfn
       character*8 aba(10000),aca(10000),baa(10000),bca(10000)
        integer iaa(10000),iba(10000),ica(10000),ida(10000),iloop(10000)
       character*512 val1(50),valp,valv
       character*15 prnvar(40)
	real invdegree
       common /basenms/adevar(200),guavar(200),cytvar(200),uravar(200)
       COMMON /HBQUA/BASE(21,20000),energy(20,20000),
     1                    yangle(20,20000)
       COMMON /HBQTYP/TYPE(20,20000),FEATURE(20,20000)
       COMMON /NAMES/NAMED(900000),RESD(900000)
       COMMON /CHAINS/CHAIND(900000)
       COMMON /ALPHA/XD(900000),YD(900000),ZD(900000)
       COMMON /BETA/IRESD(900000),ISTS(1000000),IENS(1000000),nres,
     1   allins(1000000),nforce,navar,ngvar,ncvar,nuvar,mresd(1000000)
       COMMON /DELTA/l,lresno
       COMMON /DELTAD/tmppos,tmpins
       COMMON /GAMA/ATGC(1000000)
c       common /options/cutang,cuteng,cutoff,hetatm
        common /options/cutang,cuteng,cutoff,hetatm,modelno
       common /residue/resid(900000)
       common /occur/occ(900000)
       common /num/kresd(1000000),prd(1000000)
       common /cum/pchaind(1000000),pcd(1000000),molid
       common /str/strseq(1000000)
	data ans/'No ','Yes'/
c	write(*,*) "up to this also working"
	ierror=0
        degree = 180.0/3.14159
        invdegree=3.1459/180.0
        cutoff=14.44
	hbdst=3.8
	cutang = 120.0
        cuteng = 1.80
	chohb = 1
	hetatm=0
	sugarbp = 1
	molid='    '
	nforce = 1
	oligo = 0
	nmr = 0
        cif = 0
        modelno = 1
	npoorat = 0
        nocor = 0
        nofasta = 1
        nodat = 0
        nohlx = 0     
        nonup = 0
	nodbn = 0
        nocsv = 0
	nevalue = 0
        raretauto=.false.
        if(narg.eq.0)then
	  write(6,9999)
	call showoptn(ierror)
9999	format(/'Invalid Input'/' [USAGE]:  
     1BPFIND -[options] filename'/'See BPFind_README.txt for HELP')
	ierror=1
	return
	endif
	i=0
        do while(i.lt.narg)
	  i=i+1
c	  write(*,*) "in main", val1(i)
           if(val1(i)(1:1).eq.'-') then
	    valp=val1(i)(2:20)
            if(valp.eq.'HD'.or.valp.eq.'hd'.or.valp.eq.'Hd'.or.valp.
     1       eq.'hD')then
	       i=i+1
	       valv = val1(i)
               read(valv,*,ERR=210) hbdst
	       cutoff=hbdst*hbdst
	   else if(valp.eq.'VA'.or.valp.eq.'Va'.or.valp.eq.'va'.or.valp.
     1       eq.'vA')then
	       i=i+1
	       valv = val1(i)
               read(valv,*,ERR=210) cutang
	   else if(valp.eq.'EN'.or.valp.eq.'en'.or.valp.eq.'En'.or.valp.
     1       eq.'eN')then
	       i=i+1
	       valv = val1(i)
               read(valv,*,ERR=210) cuteng
	   else if(valp.eq.'ML'.or.valp.eq.'ml'.or.valp.eq.'Ml'.or.valp.
     1       eq.'mL') then
	       i=i+1
	       valv=val1(i)
	       read(valv,*,ERR=210) molid 
           else if(valp.eq.'md'.or.valp.eq.'MD'.or.valp.eq.'Md'.or.valp.
     1        eq.'mD') then
                nmr=1
                i=i+1
                valv=val1(i)
                read(valv,*,ERR=210) modelno
	   else if(valp.eq.'HT'.or.valp.eq.'ht'.or.valp.eq.'Ht'.or.valp.
     1       eq.'hT') then
	         hetatm=1
           else if(valp.eq.'C1'.or.valp.eq.'c1') then
	         nevalue=1
	   else if(valp.eq.'CH'.or.valp.eq.'ch'.or.valp.eq.'Ch'.or.valp.
     1       eq.'cH') then
	       chohb=0
	   else if(valp.eq.'SG'.or.valp.eq.'sg'.or.valp.eq.'Sg'.or.valp.
     1       eq.'sG') then
	       sugarbp=0
	   else if(valp.eq.'AB'.or.valp.eq.'ab'.or.valp.eq.'Ab'.or.valp.
     1       eq.'aB') then
	       nforce = 0
	   else if(valp.eq.'OL'.or.valp.eq.'ol'.or.valp.eq.'Ol'.or.valp.
     1       eq.'oL') then
	       oligo=1
	    else if(valp.eq.'NMR'.or.valp.eq.'nmr'.or.valp.eq.'Nmr')then
	       nmr = 1
            else if(valp.eq.'NODAT'.or.valp.eq.'nodat'.or.valp.eq.
     1       'NoDat'.or.valp.eq.'NOdat'.or.valp.eq.'NoDAT') then
               nodat = 1
	    else if(valp.eq.'NONUP'.or.valp.eq.'nonup'.or.valp.eq.
     1       'NoNup'.or.valp.eq.'NoNUP'.or.valp.eq.'NOnup') then
               nonup = 1
	    else if(valp.eq.'NOHLX'.or.valp.eq.'nohlx'.or.valp.eq.
     1        'NoHlx'.or.valp.eq.'NoHLX'.or.valp.eq.'NOhlx') then
	       nohlx = 0     ! temporary bhatta June 2021
	    else if(valp.eq.'NOCSV'.or.valp.eq.'nocsv'.or.valp.eq.
     1       'NoCsv'.or.valp.eq.'NoCSV'.or.valp.eq.'NOcsv') then
               nocsv = 1
	    
	    else if(valp.eq.'FASTA'.or.valp.eq.'fasta'.or.valp.eq.
     1       'FAsta'.or.valp.eq.'FASta'.or.valp.eq.'FaSta') then
               nofasta = 0

  	    elseif(valp.eq.'NOCOR'.or.valp.eq.'nocor'.or.valp.eq.
     1       'NOcor'.or.valp.eq.'NOCor'.or.valp.eq.'NoCor') then
	        nocor = 1

	    else if(valp.eq.'CIF'.or.valp.eq.'cif'.or.valp.eq.'Cif'.or.
     1       valp.eq.'mmcif'.or.valp.eq.'mmCIF'.or.valp.eq.'MMCIF')then
               cif = 1
            else
               write(6,9999)
	       call showoptn(ierror)
	       return
            endif
          else 
            filenm=val1(i)
	    storedfn=filenm
          endif
        enddo
C
C Reading various types of names (three letter codes) of different modified bases.
C These files can be updated to accomodate new residue names of newly detected
C bases.
C
         call getenv('NUCLEIC_ACID_DIR',filen2)
c        filen2='../PROGRAMS'
        filen2=trim(filen2)//'/AdeVariants.name'
        open(unit=14,file=filen2,status='OLD',ERR=256)
c        open(unit=14,file='/mnt/f/Work/bpfind/AdeVariants.name')
        i=1
        do while (i.le.200)
           read(14,15,end=115) adevar(i)
           i=i+1
        enddo
115     navar=i-1
        close(unit=14)
        call getenv('NUCLEIC_ACID_DIR',filen2)
c        filen2='../PROGRAMS'
        filen2=trim(filen2)//'/GuaVariants.name'
        open(unit=14,file=filen2,status='OLD',ERR=256)
c        open(unit=14,file='/mnt/f/Work/bpfind/GuaVariants.name')
        i=1
        do while(i.le.200)
           read(14,15,end=114) guavar(i)
           i=i+1
        enddo
114     ngvar=i-1
        close(unit=14)
        call getenv('NUCLEIC_ACID_DIR',filen2)
c        filen2='../PROGRAMS'
        filen2=trim(filen2)//'/CytVariants.name'
        open(unit=14,file=filen2,status='OLD',ERR=256)
c        open(unit=14,file='/mnt/f/Work/bpfind/CytVariants.name')
        i=1
        do while(i.le.200)
           read(14,15,end=117) cytvar(i)
           i=i+1
        enddo
117     ncvar=i-1
        close(unit=14)
        call getenv('NUCLEIC_ACID_DIR',filen2)
c        filen2='../PROGRAMS'
        filen2=trim(filen2)//'/UraVariants.name'
        open(unit=14,file=filen2,status='OLD',ERR=256)
c        open(unit=14,file='/mnt/f/Work/bpfind/UraVariants.name')
        i=1
        do while(i.le.200)
           read(14,15,end=118) uravar(i)
           i=i+1
        enddo
118     continue
	close(unit=14)
        nuvar=i-1
15     format(a3)

c	write(*,*) 'Options used are:'
c	write(*,*) '        cutoff Distance^2=',cutoff
c	write(*,*) '        CutAngle=',cutang
c	write(*,*) '        CutEnergy=',cuteng
c	write(*,*) '        CH...O HB: ',ans(chohb+1)
c	write(*,*) '        Sugar HB: ',ans(sugarbp+1)
c	write(*,*) '        Choise of HETATM entries: ',ans(hetatm+1)
c	write(*,*) '        Option as Oligonucleotide: ',ans(oligo+1)
c	if(molid.ne.'    ') write(*,*) '        Mol-ID: ',molid
        if(nmr.eq.1) write(*,*) '        Analysis of first NMR model'
        if(modelno.ne.1) then
c           write(*,*) '        Analysis of',modelno,'rd model of NMR'
        endif
	kkforce = nforce + 2
	if(kkforce.gt.2) kkforce=kkforce-2
c      	write(*,*) '        Pairing restriction: ',ans(kkforce)
        title=filenm
        nnf=index(filenm,'.')
        open(unit=4,file=storedfn,status='old',err=211)
        line=filenm
	if(nocsv.eq.0) then
           line(nnf:nnf+12)='_pairing.csv'
           open(unit=53,file=line)
	   write(53,753)
	endif
753	format('#    Description of each column is the following:'/
     1 '# 1. PDB ID'/'# 2. Serial number of the cleaned residues'/
     2 '# 3. Author defined Residue number as in PDB or CIF file'/
     3 '# 4. Residue names'/'# 5. PDB_Ins_code'/'# 6. Chain name'/
     4 '# 7. Serial number of the paired residue. Blank if unpaired'/
     5 '# 8. Residue number of the paired base'/
     6 '# 9. Residue name of the paired residue'/
     6 '# 10. PDB_Ins Code of the paired residue'/
     7 '# 11. Chain name of the paired residue'/
     8 '# 12. Base pair type'/'# 13. Base pair orientation'/
     9 '# 14. These all are primery base pairs'/
     9 '# 15. Quality of base pairing as calculated by E-value'/
     1 '# 16. Residue Serial number of the base paired to the base ' 
     1 'defined in Column 2 as base multiplet. Blank if unpaired'/
     3 '# 17 to 22 are similar to columns 3 to 13'/
     4 '# 23. TP stands for Triplet Pair and BF stands for Bifurcated '
     4 'pair'/
     5 '# 24. Residue serial number of the base paired to that given in'
     5 'column 2'/
     6 '# Information of the base in Base Multiplet'/
     7 '#----------------------------------------------------')

        nmpass=filenm
        nmpass(nnf:nnf+4)='.out '
        open(unit=52,file=nmpass)
	write(52,6) filenm(1:index(filenm,'.')-1)
	if(hetatm.eq.1) write(52,375)
	if(chohb.eq.0) write(52,376)
	if(sugarbp.eq.0) write(52,377)
	if(nocor.eq.0) then
          filenm(nnf:nnf+8)='_rna.pdb '
          open(unit=9,file=filenm)
          filenm(nnf:nnf+8)='         '
	endif
	if(nonup.eq.0) then
          filenm(nnf:nnf+4)='.nup '
          open(unit=62,file=filenm)
	endif
	if(nodat.eq.0) then
          filenm(nnf:nnf+4)='.dat '
          open(unit=79,file=filenm)
	endif
	filenm(nnf:nnf+4)='.dbn'
	open(unit=80,file=filenm)
	filenm(nnf:nnf+6)='.bpseq'
	open(unit=81,file=filenm)
	if(nofasta.eq.0) then
           filenm(nnf:nnf+6)='.fasta'
           open(unit=77,file=filenm)
	endif
	if(nocsv.eq.0) then
           filenm(nnf:nnf+14)='_structure.csv'
           open(unit=78,file=filenm)
	endif
        l=1
        lresno=-9999
	tmppos=' ' 
        nres=0
        ists(1)=1
        IT = 0
        KKTI = 0
        if(cif.eq.1) then            !changes the position DM March 30,2017
	  write(52,275) ans(chohb+1)
	  write(52,276) ans(sugarbp+1)
	  write(52,277) ans(hetatm+1)
	  write(52,278) ans(oligo+1)
	  if(molid.ne.'    ') write(52,279) molid
279     format('#HEADER   Mol-ID: ',a4)
            if(nmr.eq.1) write(52,75) modelno
            write(52,73)
            write(52,74)hbdst,cutang,cuteng
            call readcif(npoorat,modelfnd)
        endif
        if(modelfnd.eq.0.or.cif.eq.0) then
        rewind(unit=4)
        do while(it.eq.0)
          READ(4,5) TITLE
            IF(TITLE(1:6).EQ.'HEADER') THEN
              KKTI = KKTI + 1
              ttl(KKTI) = TITLE
            ELSE IF(TITLE(1:6).EQ.'COMPND') THEN
              KKTI = KKTI + 1
              ttl(KKTI) = TITLE
            ELSE IF(TITLE(1:6).EQ.'AUTHOR') THEN
              KKTI = KKTI + 1
              ttl(KKTI) = TITLE
            ELSE IF(TITLE(1:4).EQ.'JRNL') THEN
              KKTI = KKTI + 1
              ttl(KKTI) = TITLE
            END IF
             IF(TITLE(1:4).EQ.'ATOM') IT = 1
             IF(KKTI.GT.49) IT = 1
	enddo
	rewind(unit=4)
	do it=1,kkti
	  write(52,255) ttl(it)
	  if(nocor.eq.0) then
	    write(9,255) ttl(it)
	  endif
	enddo
	write(52,6) filenm(1:index(filenm,'.')-1)
        write(52,73)
        write(52,74)hbdst,cutang,cuteng
C	write(52,275) ans(chohb+1)
C	write(52,276) ans(sugarbp+1)
C	write(52,277) ans(hetatm+1)
C	write(52,278) ans(oligo+1)
          if(nmr.eq.1) write(52,75) modelno
C
C  Going to read coordinates from PDB formatted file
C
        n=1
        do while(n.eq.1)
c        do 2 n=1,900000
          read(4,5,end=99)line
          if(nmr.eq.1.and.line(1:5).eq.'MODEL') read(line,79) nmrmdl
          if(nmrmdl.eq.modelno.or.nmr.eq.0) then
            do 3 kk=1,900000
              read(4,5,end=99)line
              if((line(1:4).eq.'ATOM'.or.(hetatm.eq.1.and.line(1:6).eq.
     1 'HETATM')).and.((molid.eq.' ').or.(line(22:22).eq.molid)))then
                read(line,4) dummynm,resn,chan,nores,posins,xa,ya,za,ocp
                if(dummynm(1:1).eq.' ') then
                  atmnam(1:4)=dummynm(2:5)
                else
                  atmnam(1:4)=dummynm(1:4)
                endif
                if(atmnam(4:4).eq.'A') atmnam(4:4)=' '				!what is it?
	        if(atmnam(4:4).eq.'B') write(6,*) line				!
                do j=1,2
	          if(resn(3:3).eq.' ') then
                   resn(2:3)=resn(1:2)
	           resn(1:1)=' '
	          endif
	        enddo
                chan2=chan//'   '
                call hypothesis(atmnam,resn,chan2,nores,posins,xa,
     1 ya,za,ocp,npoorat,'?')
              endif
	       if(line(1:6).eq.'ENDMDL') goto 99
3       continue
!==============================================
7          continue
	  endif
c2       continue
        enddo
99      continue
C
C   Finished READING coordinates from PDB formatted file
C
        endif
        if(npoorat.gt.1) then
          write(52,338) npoorat
        endif
        iens(nres)=l-1
338     format('#='/'#HEADER ====Number of POORLY defined and REJECTED'
     1 ' atoms:',i5,' ===='/'#=')
!-------------------------------------------------------
        do i=1,nres
            kresd(i)=mresd(i)
c           write(6,*) 'readcif :',i,mresd(i),ists(i),iens(i),kresd(i)
c           do l=ists(i),iens(i)
c            kresd(l)=mresd(l)
c              write(6,110)l,named(l),resd(l),chaind(l),kresd(i),
c     1                      xd(l),yd(l),zd(l),occ(l)
c           end do
         end do
110     format('ATOM',1X,I6,2X,A4,A3,1X,A1,I8,4X,3F8.3,2x,f4.2,F6.2)
!------------------------------------------------------
c        iens(nres)=l
        natd = l-1
        call cleanpdb(filenm,nocor,nofasta)
        write(52,160)
        write(52,77)
        write(52,60)
        write(52,61)
        write(52,62)
        write(52,63)
        write(52,163)
        write(52,64)
        write(52,65)
        write(52,66)
        write(52,67)
        write(52,164)
        write(52,68)
        write(52,69)
        write(52,70)
        write(52,71)
        write(52,72)
        do i=1,nres
         base(1,i)=0
        enddo
c        write(*,*) 'NRES=',nres
        do i = 1,nres               ! For Residue serial no. in NUPARM formalism

c*******************************ADE-ADE********************************

c       ADE-ADE (WWT) HB betn. N6A-N1A & N1A-N6A
        call findpair(i,1,'N7  ','N6  ','N1  ','N9  ',1,'N9  ','N1  ',
     1  'N6  ','N7  ','T','W:W','W:W',0)

c       N6G-ADE (WWT) HB betn. N6A-N1A & N1A-N6A
        call findpair(i,5,'N7  ','N6  ','N1  ','N9  ',1,'N9  ','N1  ',
     1  'N6  ','N7  ','T','W:W','W:W'
     2 ,0)

c       N6G-N6G (WWT) HB betn. N6A-N1A & N1A-N6A
        call findpair(i,5,'N7  ','N6  ','N1  ','N9  ',5,'N9  ','N1  ',
     1  'N6  ','N7  ','T','W:W','W:W',0)

c       ADE-ADE (WHT) HB betn. N6A-N7A & N1A-N6A
        call findpair(i,1,'N7  ','N6  ','N7  ','N9  ',1,'N9  ','N1  ',
     1  'N6  ','N1  ','T','W:H','H:W',0)

c       ADE-N6G (WHT) HB betn. N6A-N7A & N1A-N6A
        call findpair(i,1,'N7  ','N6  ','N7  ','N9  ',5,'N9  ','N1  ',
     1  'N6  ','N1  ','T','W:H','H:W',0)

c       N6G-ADE (WHT) HB betn. N6A-N7A & N1A-N6A
        call findpair(i,5,'N7  ','N6  ','N7  ','N9  ',1,'N9  ','N1  ',
     1  'N6  ','N1  ','T','W:H','H:W',0)

c       N6G-N6G (WHT) HB betn. N6A-N7A & N1A-N6A
        call findpair(i,5,'N7  ','N6  ','N7  ','N9  ',5,'N9  ','N1  ',
     1  'N6  ','N1  ','T','W:H','H:W',0)

c       ADE-ADE (HHT) HB betn. N7A-N6A & N6A-N7A
        call findpair(i,1,'N9  ','N7  ','N6  ','N1  ',1,'N1  ','N6  ',
     1  'N7  ','N9  ','T','H:H','H:H',0)

c       ADE-N6G (HHT) HB betn. N7A-N6A & N6A-N7A
        call findpair(i,1,'N9  ','N7  ','N6  ','N1  ',5,'N1  ','N6  ',
     1  'N7  ','N9  ','T','H:H','H:H',0)

c       N6G-N6G (HHT) HB betn. N7A-N6A & N6A-N7A
        call findpair(i,5,'N9  ','N7  ','N6  ','N1  ',5,'N1  ','N6  ',
     1  'N7  ','N9  ','T','H:H','H:H',0)

c*******************************ADE-GUA********************************

c       ADE-GUA (WWC) HB betn. N1A-N1G & N6A-O6G 
        call findpair(i,2,'N9  ','N1  ','N1  ','N9  ',1,'N7  ','O6  ',
     1  'N6  ','N7  ','C','W:W','W:W',0)

c       ADE-QUO (WWC) HB betn. N1A-N1G & N6A-O6G 
        call findpair(i,7,'N9  ','N1  ','N1  ','N9  ',1,'C7  ','O6  ',
     1  'N6  ','N7  ','C','W:W','W:W',0)

c       N6G-GUA (WWC) HB betn. N1A-N1G & N6A-O6G 
        call findpair(i,2,'N9  ','N1  ','N1  ','N9  ',5,'N7  ','O6  ',
     1  'N6  ','N7  ','C','W:W','W:W',0)

c       N6G-QUO (WWC) HB betn. N1A-N1G & N6A-O6G 
        call findpair(i,7,'N9  ','N1  ','N1  ','N9  ',5,'C7  ','O6  ',
     1  'N6  ','N7  ','C','W:W','W:W',0)

c       ADE-GUA (+HC) A(+N1) HB betn. N1A-N7G & N6A-O6G 
        call findpair(i,2,'N9  ','N7  ','N1  ','N9  ',1,'N1  ','O6  ',
     1  'N6  ','N7  ','C','H:+','+:H',0)

c       N6G-GUA (+HC) A(+N1) HB betn. N1A-N7G & N6A-O6G 
        call findpair(i,2,'N9  ','N7  ','N1  ','N9  ',5,'N1  ','O6  ',
     1  'N6  ','N7  ','C','H:+','+:H',0)

c       ADE-GUA (+HT) A(N1+) HB betn. N6A-N7G & N1A-O6G 
        call findpair(i,2,'N9  ','N7  ','N6  ','N7  ',1,'N1  ','O6  ',
     1  'N1  ','N9  ','T','H:+','+:H',0)

c       N6G-GUA (+HT) A(N1+) HB betn. N6A-N7G & N1A-O6G 
        call findpair(i,2,'N9  ','N7  ','N6  ','N7  ',5,'N1  ','O6  ',
     1  'N1  ','N9  ','T','H:+','+:H',0)

c       ADE-GUA (WST) HB betn. N1A-N2G & N6A-N3G
        call findpair(i,1,'N9  ','N1  ','N2  ','N1  ',2,'N7  ','N6  ',
     1  'N3  ','C6  ','T','W:S','S:W',0)

c       ADE-QUO (WST) HB betn. N1A-N2G & N6A-N3G
        call findpair(i,1,'N9  ','N1  ','N2  ','N1  ',7,'N7  ','N6  ',
     1  'N3  ','C6  ','T','W:S','S:W',0)

c       N6G-GUA (WST) HB betn. N1A-N2G & N6A-N3G
        call findpair(i,5,'N9  ','N1  ','N2  ','N1  ',2,'N7  ','N6  ',
     1  'N3  ','C6  ','T','W:S','S:W',0)

c       N6G-QUO (WST) HB betn. N1A-N2G & N6A-N3G
        call findpair(i,5,'N9  ','N1  ','N2  ','N1  ',7,'N7  ','N6  ',
     1  'N3  ','C6  ','T','W:S','S:W',0)

c       ADE-N6G (WST) HB betn. N1A-N2G & N6A-N3G
        call findpair(i,1,'N9  ','N1  ','N2  ','N1  ',5,'N7  ','N6  ',
     1  'N3  ','C6  ','T','W:S','S:W',0)

c       N6G-N6G (WST) HB betn. N1A-N2G & N6A-N3G
        call findpair(i,5,'N9  ','N1  ','N2  ','N1  ',5,'N7  ','N6  ',
     1  'N3  ','C6  ','T','W:S','S:W',0)

c       ADE-GUA (HWC) HB betn. N6A-O6G & N7A-N1G
        call findpair(i,2,'N7  ','O6  ','N6  ','N1  ',1,'N9  ','N1  ',
     1  'N7  ','C4  ','C','W:H','H:W',0)

c       ADE-QUO (HWC) HB betn. N6A-O6G & N7A-N1G
        call findpair(i,7,'C7  ','O6  ','N6  ','N1  ',1,'N9  ','N1  ',
     1  'N7  ','C4  ','C','W:H','H:W',0)

c       N6G-GUA (HWC) HB betn. N6A-O6G & N7A-N1G
        call findpair(i,2,'N7  ','O6  ','N6  ','N1  ',5,'N9  ','N1  ',
     1  'N7  ','C4  ','C','W:H','H:W',0)

c       N6G-QUO (HWC) HB betn. N6A-O6G & N7A-N1G
        call findpair(i,7,'C7  ','O6  ','N6  ','N1  ',5,'N9  ','N1  ',
     1  'N7  ','C4  ','C','W:H','H:W',0)

c       ADE-GUA (HST) HB betn. N6A-N3G & N7A-N2G 
        call findpair(i,2,'C6  ','N3  ','N6  ','N1  ',1,'N1  ','N2  ',
     1  'N7  ','N9  ','T','S:H','H:S',0)

c       ADE-QUO (HST) HB betn. N6A-N3G & N7A-N2G 
        call findpair(i,7,'C6  ','N3  ','N6  ','N1  ',1,'N1  ','N2  ',
     1  'N7  ','N9  ','T','S:H','H:S',0)

c       ADE-N6G (HST) HB betn. N6A-N3G & N7A-N2G 
        call findpair(i,5,'C6  ','N3  ','N6  ','N1  ',1,'N1  ','N2  ',
     1  'N7  ','N9  ','T','S:H','H:S',0)

c       N6G-GUA (HST) HB betn. N6A-N3G & N7A-N2G 
        call findpair(i,2,'C6  ','N3  ','N6  ','N1  ',5,'N1  ','N2  ',
     1  'N7  ','N9  ','T','S:H','H:S',0)

c       N6G-QUO (HST) HB betn. N6A-N3G & N7A-N2G 
        call findpair(i,7,'C6  ','N3  ','N6  ','N1  ',5,'N1  ','N2  ',
     1  'N7  ','N9  ','T','S:H','H:S',0)

c       N6G-N6G(HST) HB betn. N6A-N3G & N7A-N2G 
        call findpair(i,5,'C6  ','N3  ','N6  ','N1  ',5,'N1  ','N2  ',
     1  'N7  ','N9  ','T','S:H','H:S',0)

c*******************************ADE-CYT********************************

c       ADE-CYT (+WC) A(+N1) HB betn. N1A-O2C & N6A-N3C
        call findpair(i,1,'N9  ','N1  ','O2  ','C1* ',3,'N7  ','N6  ',
     1  'N3  ','C6  ','C','+:W','W:+',0)

c       N6G-CYT (+WC) A(+N1) HB betn. N1A-O2C & N6A-N3C
        call findpair(i,5,'N9  ','N1  ','O2  ','C1* ',3,'N7  ','N6  ',
     1  'N3  ','C6  ','C','+:W','W:+',0)

c       ADE-CYT (WWT) HB betn. N1A-N4C & N6A-N3C
        call findpair(i,1,'N9  ','N1  ','N4  ','C5  ',3,'N7  ','N6  ',
     1  'N3  ','C6  ','T','W:W','W:W',0)

c       N6G-CYT (WWT) HB betn. N1A-N4C & N6A-N3C
        call findpair(i,5,'N9  ','N1  ','N4  ','C5  ',3,'N7  ','N6  ',
     1  'N3  ','C6  ','T','W:W','W:W',0)

c       ADE-CYT (W+T) (C+&AW) HB betn. N1A-N3C & N6A-O2C
        call findpair(i,3,'C5  ','N3  ','N1  ','C4  ',1,'N1  ','O2  ',
     1  'N6  ','N7  ','T','+:W','W:+',0)

c       N6G-CYT (W+T) (C+&AW) HB betn. N1A-N3C & N6A-O2C
        call findpair(i,3,'C5  ','N3  ','N1  ','C4  ',5,'N1  ','O2  ',
     1  'N6  ','N7  ','T','+:W','W:+',0)

c       ADE-CYT (HWT) HB betn. N7A-N4C & N6A-N3C 
        call findpair(i,3,'C5  ','N4  ','N7  ','N9  ',1,'C6  ','N3  ',
     1  'N6  ','N1  ','T','W:H','H:W',0)

c       N6G-CYT (HWT) HB betn. N7A-N4C & N6A-N3C 
        call findpair(i,3,'C5  ','N4  ','N7  ','N9  ',5,'C6  ','N3  ',
     1  'N6  ','N1  ','T','W:H','H:W',0)

c       ADE-CYT (H+T) HB betn. N7A-N3C & N6A-O2C
c        call findpair(i,3,'C6  ','N3  ','N7  ','N9  ',1,'C1* ','O2  ',
c     1  'N6  ','N1  ','T','+:H','H:+',0)    Very similar to HWT, Feb. 22, 2010

c       N6G-CYT (H+T) HB betn. N7A-N3C & N6A-O2C
        call findpair(i,3,'C6  ','N3  ','N7  ','N9  ',5,'C1* ','O2  ',
     1  'N6  ','N1  ','T','+:H','H:+',0)

c*******************************ADE-URA********************************

c       ADE-URA (WWC) Canonical HB betn. N1A-N3U & N6A-O4U
        call findpair(i,1,'N9  ','N1  ','N3  ','C6  ',4,'N7  ','N6  ',
     1  'O4  ','C5  ','C','W:W','W:W',0)

c       N6G-URA (WWC) Canonical HB betn. N1A-N3U & N6A-O4U
        call findpair(i,5,'N9  ','N1  ','N3  ','C6  ',4,'N7  ','N6  ',
     1  'O4  ','C5  ','C','W:W','W:W',0)

c       ADE-PSU (WWC) HB betn. N1A-N3U & N6A-O2U  
        call findpair(i,1,'N9  ','N1  ','N3  ','C6  ',6,'N7  ','N6  ',
     1  'O2  ','N1  ','C','W:W','W:W',0)

c       N6G-PSU (WWC) HB betn. N1A-N3U & N6A-O2U  
        call findpair(i,5,'N9  ','N1  ','N3  ','C6  ',6,'N7  ','N6  ',
     1  'O2  ','N1  ','C','W:W','W:W',0)

c       ADE-URA (WWT) HB betn. N1A-N3U & N6A-O2U  
        call findpair(i,1,'N9  ','N1  ','N3  ','C6  ',4,'N7  ','N6  ',
     1  'O2  ','C1* ','T','W:W','W:W',0)

c       N6G-URA (WWT) HB betn. N1A-N3U & N6A-O2U  
        call findpair(i,5,'N9  ','N1  ','N3  ','C6  ',4,'N7  ','N6  ',
     1  'O2  ','C1* ','T','W:W','W:W',0)

c       ADE-PSU (WWT) HB betn. N1A-N3U & N6A-O4U
        call findpair(i,1,'N9  ','N1  ','N3  ','C6  ',6,'N7  ','N6  ',
     1  'O4  ','C1* ','T','W:W','W:W',0)

c       N6G-PSU (WWT) HB betn. N1A-N3U & N6A-O4U
        call findpair(i,5,'N9  ','N1  ','N3  ','C6  ',6,'N7  ','N6  ',
     1  'O4  ','C1* ','T','W:W','W:W',0)

c       ADE-PSU (WHC) HB betn. N1A-N1U & N6A-O2U
        call findpair(i,1,'N9  ','N1  ','N1  ','C4  ',6,'N7  ','N6  ',
     1  'O2  ','N3  ','C','W:H','H:W',0)

c       N6G-PSU (WHC) HB betn. N1A-N1U & N6A-O2U
        call findpair(i,5,'N9  ','N1  ','N1  ','C4  ',6,'N7  ','N6  ',
     1  'O2  ','N3  ','C','W:H','H:W',0)

c       ADE-URA (HWC) HB betn. N6A-O4U & N7A-N3U 
        call findpair(i,1,'N1  ','N6  ','O4  ','C5  ',4,'N9  ','N7  ',
     1  'N3  ','C6  ','C','H:W','W:H',0)

c       N6G-URA (HWC) HB betn. N6A-O4U & N7A-N3U 
        call findpair(i,5,'N1  ','N6  ','O4  ','C5  ',4,'N9  ','N7  ',
     1  'N3  ','C6  ','C','H:W','W:H',0)

c       ADE-PSU (HWC) HB betn. N7A-N3U & N6A-O2U 
        call findpair(i,1,'N9  ','N7  ','N3  ','C6  ',6,'N1  ','N6  ',
     1  'O2  ','N1  ','C','H:W','W:H',0)

c       N6G-PSU (HWC) HB betn. N7A-N3U & N6A-O2U 
        call findpair(i,5,'N9  ','N7  ','N3  ','C6  ',6,'N1  ','N6  ',
     1  'O2  ','N1  ','C','H:W','W:H',0)

c       ADE-URA (HWT) HB betn. N7A-N3U & N6A-O2U 
        call findpair(i,1,'N9  ','N7  ','N3  ','C6  ',4,'N1  ','N6  ',
     1  'O2  ','C1* ','T','H:W','W:H',0)

c       N6G-URA (HWT) HB betn. N7A-N3U & N6A-O2U 
        call findpair(i,5,'N9  ','N7  ','N3  ','C6  ',4,'N1  ','N6  ',
     1  'O2  ','C1* ','T','H:W','W:H',0)

c       ADE-PSU (HWT) HB betn. N6A-O4U & N7A-N3U 
        call findpair(i,1,'N1  ','N6  ','O4  ','C1* ',6,'N9  ','N7  ',
     1  'N3  ','C6  ','T','H:W','W:H',0)

c       N6G-PSU (HWT) HB betn. N6A-O4U & N7A-N3U 
        call findpair(i,5,'N1  ','N6  ','O4  ','C1* ',6,'N9  ','N7  ',
     1  'N3  ','C6  ','T','H:W','W:H',0)

c       ADE-PSU (HHC) HB betn. N7A-N1U & N6A-O2U
        call findpair(i,1,'C4  ','N7  ','N1  ','C4  ',6,'N1  ','N6  ',
     1  'O2  ','N3  ','C','H:H','H:H',0)

c       N6G-PSU (HHC) HB betn. N7A-N1U & N6A-O2U
        call findpair(i,5,'C4  ','N7  ','N1  ','C4  ',6,'N1  ','N6  ',
     1  'O2  ','N3  ','C','H:H','H:H',0)

c*******************************GUA-GUA********************************

c       GUA-GUA (WWT) HB betn. O6G-N1G & N1G-O6G 
        call findpair(i,2,'N7  ','O6  ','N1  ','N9  ',2,'N9  ','N1  ',
     1  'O6  ','N7  ','T','W:W','W:W',0)

c       GUA-QUO (WWT) HB betn. O6G-N1G & N1G-O6G 
        call findpair(i,2,'N7  ','O6  ','N1  ','N9  ',7,'N9  ','N1  ',
     1  'O6  ','C7  ','T','W:W','W:W',0)

c       QUO-QUO (WWT) HB betn. O6G-N1G & N1G-O6G 
        call findpair(i,7,'C7  ','O6  ','N1  ','N9  ',7,'N9  ','N1  ',
     1  'O6  ','C7  ','T','W:W','W:W',0)

c       GUA-GUA (WHC) HB betn. N1G-O6G & N2G-N7G
        call findpair(i,2,'N9  ','N1  ','O6  ','N1  ',2,'N3  ','N2  ',
     1  'N7  ','N9  ','C','W:H','H:W',0)

c       QUO-GUA (WHC) HB betn. N1G-O6G & N2G-N7G
        call findpair(i,7,'N9  ','N1  ','O6  ','N1  ',2,'N3  ','N2  ',
     1  'N7  ','N9  ','C','W:H','H:W',0)

c       N6G-GUA (WHC) HB betn. N1G-O6G & N2G-N7G
        call findpair(i,5,'N9  ','N1  ','O6  ','N1  ',2,'N3  ','N2  ',
     1  'N7  ','N9  ','C','W:H','H:W',0)

c       GUA-GUA (WHT) HB betn. N1G-N7G & N2G-O6G 
        call findpair(i,2,'N9  ','N7  ','N1  ','N9  ',2,'N1  ','O6  ',
     1  'N2  ','N3  ','T','H:W','W:H',0)

c       QUO-GUA (WHT) HB betn. N1G-N7G & N2G-O6G 
        call findpair(i,2,'N9  ','N7  ','N1  ','N9  ',7,'N1  ','O6  ',
     1  'N2  ','N3  ','T','H:W','W:H',0)

c       N6G-GUA (WHT) HB betn. N1G-N7G & N2G-O6G 
        call findpair(i,2,'N9  ','N7  ','N1  ','N9  ',5,'N1  ','O6  ',
     1  'N2  ','N3  ','T','H:W','W:H',0)

c       GUA-GUA (WSC) HB betn. N1G-N3G & O6G-N2G 
        call findpair(i,2,'C6  ','N3  ','N1  ','N9  ',2,'N1  ','N2  ',
     1  'O6  ','N7  ','C','S:W','W:S',0)

c       GUA-QUO (WSC) HB betn. N1G-N3G & O6G-N2G 
        call findpair(i,7,'C6  ','N3  ','N1  ','N9  ',2,'N1  ','N2  ',
     1  'O6  ','N7  ','C','S:W','W:S',0)

c       QUO-GUA (WSC) HB betn. N1G-N3G & O6G-N2G 
        call findpair(i,2,'C6  ','N3  ','N1  ','N9  ',7,'N1  ','N2  ',
     1  'O6  ','C7  ','C','S:W','W:S',0)

c       QUO-QUO (WSC) HB betn. N1G-N3G & O6G-N2G 
        call findpair(i,7,'C6  ','N3  ','N1  ','N9  ',7,'N1  ','N2  ',
     1  'O6  ','C7  ','C','S:W','W:S',0)

c       GUA-N6G (WSC) HB betn. N1G-N3G & O6G-N2G 
        call findpair(i,5,'C6  ','N3  ','N1  ','N9  ',2,'N1  ','N2  ',
     1  'O6  ','N7  ','C','S:W','W:S',0)

c       QUO-N6G (WSC) HB betn. N1G-N3G & O6G-N2G 
        call findpair(i,5,'C6  ','N3  ','N1  ','N9  ',7,'N1  ','N2  ',
     1  'O6  ','C7  ','C','S:W','W:S',0)

c       GUA-GUA (HzT) HB betn. O6G-N3G & N7G-N2G  
        call findpair(i,2,'C6  ','N3  ','O6  ','N1  ',2,'N1  ','N2  ',
     1  'N7  ','N9  ','T','z:H','H:z',0)

c       GUA-QUO (HzT) HB betn. O6G-N3G & N7G-N2G  
        call findpair(i,7,'C6  ','N3  ','O6  ','N1  ',2,'N1  ','N2  ',
     1  'N7  ','N9  ','T','z:H','H:z',0)

c       GUA-N6G (HzT) HB betn. O6G-N3G & N7G-N2G  
        call findpair(i,5,'C6  ','N3  ','O6  ','N1  ',2,'N1  ','N2  ',
     1  'N7  ','N9  ','T','z:H','H:z',0)

c       GUA-GUA (SST) HB betn. N2G-N3G & N3G-N2G
        call findpair(i,2,'N1  ','N2  ','N3  ','C6  ',2,'C6  ','N3  ',
     1  'N2  ','N1  ','T','S:S','S:S',0)

c       GUA-QUO (SST) HB betn. N2G-N3G & N3G-N2G
        call findpair(i,2,'N1  ','N2  ','N3  ','C6  ',7,'C6  ','N3  ',
     1  'N2  ','N1  ','T','S:S','S:S',0)

c       QUO-QUO (SST) HB betn. N2G-N3G & N3G-N2G
        call findpair(i,7,'N1  ','N2  ','N3  ','C6  ',7,'C6  ','N3  ',
     1  'N2  ','N1  ','T','S:S','S:S',0)

c       GUA-N6G (SST) HB betn. N2G-N3G & N3G-N2G
        call findpair(i,2,'N1  ','N2  ','N3  ','C6  ',5,'C6  ','N3  ',
     1  'N2  ','N1  ','T','S:S','S:S',0)

c       QUO-N6G (SST) HB betn. N2G-N3G & N3G-N2G
        call findpair(i,7,'N1  ','N2  ','N3  ','C6  ',5,'C6  ','N3  ',
     1  'N2  ','N1  ','T','S:S','S:S',0)

c       N6G-N6G (SST) HB betn. N2G-N3G & N3G-N2G
        call findpair(i,5,'N1  ','N2  ','N3  ','C6  ',5,'C6  ','N3  ',
     1  'N2  ','N1  ','T','S:S','S:S',0)

c*******************************GUA-CYT********************************

c       GUA-CYT (WWC) Canonical HB betn. N1G-N3C, O6G-N4C & N2G-O2C
        call findhb3(i,3,2,'C6  ','N3  ','N1  ','N9  ','C5  ','N4  ',
     1  'O6  ','N7  ','N1  ','O2  ','N2  ','N3  ','C','W:W','W:W')

c       QUO-CYT (WWC) Canonical HB betn. N1G-N3C, O6G-N4C & N2G-O2C
        call findhb3(i,3,7,'C6  ','N3  ','N1  ','N9  ','C5  ','N4  ',
     1  'O6  ','C7  ','N1  ','O2  ','N2  ','N3  ','C','W:W','W:W')

c       N6G-CYT (WWC) HB betn. N1G-N3C & N2G-O2C
        call findpair(i,3,'C6  ','N3  ','N1  ','N9  ',5,'N1  ','O2  ',
     1  'N2  ','N3  ','C','W:W','W:W',0)

c       GUA-CYT (W+C) HB betn. O6G-N3C & N1G-O2C
        call findpair(i,2,'N7  ','O6  ','N3  ','C6  ',3,'C4  ','N1  ',
     1  'O2  ','N1  ','C','W:+','+:W',0)

c       QUO-CYT (W+C) HB betn. O6G-N3C & N1G-O2C
        call findpair(i,7,'C7  ','O6  ','N3  ','C6  ',3,'C4  ','N1  ',
     1  'O2  ','N1  ','C','W:+','+:W',0)

c       GUA-CYT (WWT) HB betn. N1G-O2C & N2G-N3C
        call findpair(i,2,'N9  ','N1  ','O2  ','N1  ',3,'N3  ','N2  ',
     1  'N3  ','C6  ','T','W:W','W:W',0)

c       QUO-CYT (WWT) HB betn. N1G-O2C & N2G-N3C
        call findpair(i,7,'N9  ','N1  ','O2  ','N1  ',3,'N3  ','N2  ',
     1  'N3  ','C6  ','T','W:W','W:W',0)

c       N6G-CYT (WWT) HB betn. N1G-O2C & N2G-N3C
        call findpair(i,5,'N9  ','N1  ','O2  ','N1  ',3,'N3  ','N2  ',
     1  'N3  ','C6  ','T','W:W','W:W',0)

c       GUA-CYT (H+C) C(N3+) HB betn. O6G-N4C & N7G-N3C 
        call findpair(i,3,'C5  ','N4  ','O6  ','N1  ',2,'C6  ','N3  ',
     1  'N7  ','N9  ','C','+:H','H:+',0)

c       GUA-CYT (zWC) G(N3+) HB betn. N3G-O2C & N2G-N3C
        call findpair(i,2,'C6  ','N3  ','O2  ','C1* ',3,'N1  ','N2  ',
     1  'N3  ','C6  ','C','z:W','W:z',0)

c       QUO-CYT (zWC) G(N3+) HB betn. N3G-O2C & N2G-N3C
        call findpair(i,7,'C6  ','N3  ','O2  ','C1* ',3,'N1  ','N2  ',
     1  'N3  ','C6  ','C','z:W','W:z',0)

c       N6G-CYT (zWC) G(N3+) HB betn. N3G-O2C & N2G-N3C
        call findpair(i,5,'C6  ','N3  ','O2  ','C1* ',3,'N1  ','N2  ',
     1  'N3  ','C6  ','C','z:W','W:z',0)

c       GUA-CYT (H+T) C(N3+) HB betn. N7G-N4C & O6G-N3C 
        call findpair(i,3,'C5  ','N4  ','N7  ','N9  ',2,'C6  ','N3  ',
     1  'O6  ','N1  ','T','+:H','H:+',0)

c       GUA-CYT (SWT) HB betn. N2G-N3C & N3G-N4C 
        call findpair(i,2,'N1  ','N2  ','N3  ','C6  ',3,'C6  ','N3  ',
     1  'N4  ','C5  ','T','S:W','W:S',0)

c       QUO-CYT (SWT) HB betn. N2G-N3C & N3G-N4C 
        call findpair(i,7,'N1  ','N2  ','N3  ','C6  ',3,'C6  ','N3  ',
     1  'N4  ','C5  ','T','S:W','W:S',0)

c       N6G-CYT (SWT) HB betn. N2G-N3C & N3G-N4C 
        call findpair(i,5,'N1  ','N2  ','N3  ','C6  ',3,'C6  ','N3  ',
     1  'N4  ','C5  ','T','S:W','W:S',0)

c       GUA-CYT (S+T) C(N3+) HB betn. N3G-N3C & N2G-O2C
        call findpair(i,2,'C6  ','N3  ','N3  ','C6  ',3,'N1  ','N2  ',
     1  'O2  ','C1* ','T','S:+','+:S',0)

c       QUO-CYT (S+T) C(N3+) HB betn. N3G-N3C & N2G-O2C
        call findpair(i,7,'C6  ','N3  ','N3  ','C6  ',3,'N1  ','N2  ',
     1  'O2  ','C1* ','T','S:+','+:S',0)

c       N6G-CYT (S+T) C(N3+) HB betn. N3G-N3C & N2G-O2C
        call findpair(i,5,'C6  ','N3  ','N3  ','C6  ',3,'N1  ','N2  ',
     1  'O2  ','C1* ','T','S:+','+:S',0)

c*******************************GUA-URA********************************

c       GUA-URA (WWC) HB betn. O6G-N3U & N1G-O2U 
        call findpair(i,2,'N7  ','O6  ','N3  ','C6  ',4,'N9  ','N1  ',
     1  'O2  ','C1* ','C','W:W','W:W',0)

c       QUO-URA (WWC) HB betn. O6G-N3U & N1G-O2U 
        call findpair(i,7,'C7  ','O6  ','N3  ','C6  ',4,'N9  ','N1  ',
     1  'O2  ','C1* ','C','W:W','W:W',0)

c       GUA-PSU (WWC) HB betn. N1G-O4U & O6G-N3U 
        call findpair(i,2,'N9  ','N1  ','O4  ','C1* ',6,'N7  ','O6  ',
     1  'N3  ','C6  ','C','W:W','W:W',0)

c       QUO-PSU (WWC) HB betn. N1G-O4U & O6G-N3U 
        call findpair(i,7,'N9  ','N1  ','O4  ','C1* ',6,'C7  ','O6  ',
     1  'N3  ','C6  ','C','W:W','W:W',0)

c       GUA-URA (WWT) HB betn. N1G-O4U & O6G-N3U 
        call findpair(i,2,'N9  ','N1  ','O4  ','C5  ',4,'N7  ','O6  ',
     1  'N3  ','C6  ','T','W:W','W:W',0)

c       QUO-URA (WWT) HB betn. N1G-O4U & O6G-N3U 
        call findpair(i,7,'N9  ','N1  ','O4  ','C5  ',4,'C7  ','O6  ',
     1  'N3  ','C6  ','T','W:W','W:W',0)

c       GUA-PSU (WWT) HB betn. O6G-N3U & N1G-O2U 
        call findpair(i,2,'N7  ','O6  ','N3  ','C6  ',6,'N9  ','N1  ',
     1  'O2  ','N1  ','T','W:W','W:W',0)

c       QUO-PSU (WWT) HB betn. O6G-N3U & N1G-O2U 
        call findpair(i,7,'C7  ','O6  ','N3  ','C6  ',6,'N9  ','N1  ',
     1  'O2  ','N1  ','T','W:W','W:W',0)

c       GUA-PSU (WHC) HB betn. O6G-N1U & N1G-O2U 
        call findpair(i,2,'N9  ','O6  ','N1  ','C4  ',6,'N3  ','N1  ',
     1  'O2  ','N3  ','C','W:H','H:W',0)

c       QUO-PSU (WHC) HB betn. O6G-N1U & N1G-O2U 
        call findpair(i,7,'N9  ','O6  ','N1  ','C4  ',6,'N3  ','N1  ',
     1  'O2  ','N3  ','C','W:H','H:W',0)

c       GUA-URA (SWC) HB betn. N3G-N3U & N2G-O4U
        call findpair(i,2,'C6  ','N3  ','N3  ','C6  ',4,'N1  ','N2  ',
     1  'O4  ','C5  ','C','S:W','W:S',0)

c       QUO-URA (SWC) HB betn. N3G-N3U & N2G-O4U
        call findpair(i,7,'C6  ','N3  ','N3  ','C6  ',4,'N1  ','N2  ',
     1  'O4  ','C5  ','C','S:W','W:S',0)

c       N6G-URA (SWC) HB betn. N3G-N3U & N2G-O4U
        call findpair(i,5,'C6  ','N3  ','N3  ','C6  ',4,'N1  ','N2  ',
     1  'O4  ','C5  ','C','S:W','W:S',0)

c       GUA-PSU (SWC) HB betn. N3G-N3U & N2G-O2U  
        call findpair(i,2,'C6  ','N3  ','N3  ','C6  ',6,'N1  ','N2  ',
     1  'O2  ','N1  ','C','S:W','W:S',0)

c       QUO-PSU (SWC) HB betn. N3G-N3U & N2G-O2U  
        call findpair(i,7,'C6  ','N3  ','N3  ','C6  ',6,'N1  ','N2  ',
     1  'O2  ','N1  ','C','S:W','W:S',0)

c       N6G-PSU (SWC) HB betn. N3G-N3U & N2G-O2U  
        call findpair(i,5,'C6  ','N3  ','N3  ','C6  ',6,'N1  ','N2  ',
     1  'O2  ','N1  ','C','S:W','W:S',0)

c       GUA-URA (SWT) HB betn. N3G-N3U & N2G-O2U  
        call findpair(i,2,'C6  ','N3  ','N3  ','C6  ',4,'N1  ','N2  ',
     1  'O2  ','C1* ','T','S:W','W:S',0)

c       QUO-URA (SWT) HB betn. N3G-N3U & N2G-O2U  
        call findpair(i,7,'C6  ','N3  ','N3  ','C6  ',4,'N1  ','N2  ',
     1  'O2  ','C1* ','T','S:W','W:S',0)

c       N6G-URA (SWT) HB betn. N3G-N3U & N2G-O2U  
        call findpair(i,5,'C6  ','N3  ','N3  ','C6  ',4,'N1  ','N2  ',
     1  'O2  ','C1* ','T','S:W','W:S',0)

c       GUA-PSU (SWT) HB betn. N3G-N3U & N2G-O4U
        call findpair(i,2,'C6  ','N3  ','N3  ','C6  ',6,'N1  ','N2  ',
     1  'O4  ','C1* ','T','S:W','W:S',0)

c       QUO-PSU (SWT) HB betn. N3G-N3U & N2G-O4U
        call findpair(i,7,'C6  ','N3  ','N3  ','C6  ',6,'N1  ','N2  ',
     1  'O4  ','C1* ','T','S:W','W:S',0)

c       N6G-PSU (SWT) HB betn. N3G-N3U & N2G-O4U
        call findpair(i,5,'C6  ','N3  ','N3  ','C6  ',6,'N1  ','N2  ',
     1  'O4  ','C1* ','T','S:W','W:S',0)

c*******************************CYT-CYT********************************

c       CYT-CYT (W+C) C(+N3) HB betn. O2C-N3C & N3C-N4C
        call findpair(i,3,'C1* ','O2  ','N3  ','C6  ',3,'C6  ','N3  ',
     1  'N4  ','C5  ','C','W:+','+:W',0)

c       CYT-CYT (WWT) HB Betn. N4C-N3C & N3C-N4C
        call findpair(i,3,'C5  ','N4  ','N3  ','C6  ',3,'C6  ','N3  ',
     1  'N4  ','C5  ','T','W:W','W:W',0)

c       CYT-CYT (W+T) HB betn. N3C-N3C; N4C-O2C & O2C-N4C 
        call findhb3(i,3,3,'C6  ','N3  ','N3  ','C6  ','C5  ','N4  ',
     1  'O2  ','C1* ','C1* ','O2  ','N4  ','C5  ','T','W:+','+:W')

c*******************************CYT-URA********************************

c       CYT-URA (WWC) HB betn. N3C-N3U & N4C-O4U 
        call findpair(i,4,'C6  ','N3  ','N3  ','C6  ',3,'C5  ','O4  ',
     1  'N4  ','C5  ','C','W:W','W:W',0)

c       CYT-URA (+WC) HB betn. N3C-O4U & O2C-N3U 
        call findpair(i,3,'C6  ','N3  ','O4  ','C5  ',4,'C1* ','O2  ',
     1  'N3  ','C6  ','C','+:W','W:+',0)

c       CYT-PSU (+WC) HB betn. N3C-O2U & O2C-N3U 
        call findpair(i,6,'N1  ','O2  ','N3  ','C6  ',3,'C6  ','N3  ',
     1  'O2  ','C1* ','C','W:+','+:W',0)

c       CYT-URA (WWT) HB betn. N3C-N3U & N4C-O2U
        call findpair(i,3,'C6  ','N3  ','N3  ','C6  ',4,'C5  ','N4  ',
     1  'O2  ','C1* ','T','W:W','W:W',0)

c       CYT-PSU (WWT) HB betn. N3C-N3U & N4C-O4U 
        call findpair(i,6,'C6  ','N3  ','N3  ','C6  ',3,'C1* ','O4  ',
     1  'N4  ','C5  ','T','W:W','W:W',0)

c       CYT-URA (+WT) HB betn. N3C-O2U & O2C-N3U 
        call findpair(i,4,'C1* ','O2  ','N3  ','C6  ',3,'C6  ','N3  ',
     1  'O2  ','C1* ','T','W:+','+:W',0)

c       CYT-PSU (+WT) HB betn. N3C-O4U & O2C-N3U 
        call findpair(i,3,'C6  ','N3  ','O4  ','C1* ',6,'C1* ','O2  ',
     1  'N3  ','C6  ','T','+:W','W:+',0)

c       CYT-PSU (WHC) HB betn. N4C-O2U & O2C-N1U
        call findpair(i,3,'C5  ','N4  ','O2  ','N3  ',6,'C1* ','O2  ',
     1  'N1  ','C4  ','C','W:H','H:W',0)

c*******************************URA-URA********************************

c       URA-URA (WWC) HB betn. N3U-O4U & O2U-N3U
        call findpair(i,4,'C6  ','N3  ','O4  ','C5  ',4,'C1* ','O2  ',
     1  'N3  ','C6  ','C','W:W','W:W',0)

c       URA-PSU (WWC) 1st HB betn. O4U-N3U & N3U-O4U
        call findpair(i,4,'C5  ','O4  ','N3  ','C6  ',6,'C6  ','N3  ',
     1  'O4  ','C1* ','C','W:W','W:W',0)

c       URA-PSU (WWC) 2nd HB betn. O2U-N3U & N3U-O2U
        call findpair(i,4,'C1* ','O2  ','N3  ','C6  ',6,'C6  ','N3  ',
     1  'O2  ','N1  ','C','W:W','W:W',0)

c       PSU-PSU (WWC) HB betn. N3U-O4U & O2U-N3U
        call findpair(i,6,'C6  ','N3  ','O4  ','C1* ',6,'N1  ','O2  ',
     1  'N3  ','C6  ','C','W:W','W:W',0)

c       URA-URA (WWT) 1st HB betn. O4U-N3U & N3U-O4U
        call findpair(i,4,'C5  ','O4  ','N3  ','C6  ',4,'C6  ','N3  ',
     1  'O4  ','C5  ','T','W:W','W:W',0)

c       URA-URA (WWT) 2nd HB betn. O2U-N3U & N3U-O2U
        call findpair(i,4,'C1* ','O2  ','N3  ','C6  ',4,'C6  ','N3  ',
     1  'O2  ','C1* ','T','W:W','W:W',0)

c       URA-PSU (WWT) HB betn. N3U-O4U & O2U-N3U
        call findpair(i,4,'C6  ','N3  ','O4  ','C1* ',6,'C1* ','O2  ',
     1  'N3  ','C6  ','T','W:W','W:W',0)

c       PSU-PSU (WWT) 1st HB betn. O4U-N3U & N3U-O4U
        call findpair(i,6,'C1* ','O4  ','N3  ','C6  ',6,'C6  ','N3  ',
     1  'O4  ','C1* ','T','W:W','W:W',0)

c       PSU-PSU (WWT) 2nd HB betn. O2U-N3U & N3U-O2U
        call findpair(i,6,'N1  ','O2  ','N3  ','C6  ',6,'C6  ','N3  ',
     1  'O2  ','N1  ','T','W:W','W:W',0)

c       URA-PSU (WHC) HB betn. O2U-N1U & N3U-O2U
        call findpair(i,4,'N1  ','O2  ','N1  ','C5  ',6,'C5  ','N3  ',
     1  'O2  ','N3  ','C','W:H','H:W',0)

c       PSU-PSU (WHC) HB betn. O4U-N1U & N3U-O2U
        call findpair(i,6,'C1* ','O4  ','N1  ','C4  ',6,'N1  ','N3  ',
     1  'O2  ','N3  ','C','W:H','H:W',0)

c       URA-PSU (WHT) HB betn. O4U-N1U & N3U-O2U
        call findpair(i,4,'C5  ','O4  ','N1  ','C4  ',6,'N1  ','N3  ',
     1  'O2  ','N3  ','T','W:H','H:W',0)

c       PSU-PSU (WHT) HB betn. O2U-N1U & N3U-O2U
        call findpair(i,6,'N1  ','O2  ','N1  ','C5  ',6,'C5  ','N3  ',
     1  'O2  ','N3  ','T','W:H','H:W',0)

c**********************************************************************

c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c	%%%%%           C-H...O mediated interactions           %%%%%
c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if(chohb.eq.1) then

c*******************************ADE-ADE********************************

c       ADE-ADE (wwC) HB betn. N1A-N6A & C2A-N1A
        call findpair(i,1,'C5  ','N1  ','N6  ','N7  ',1,'C4  ','C2  ',
     1  'N1  ','C4  ','C','w:w','w:w',0)

c       ADE-N6G (wwC) HB betn. N1A-N6A & C2A-N1A
        call findpair(i,1,'C5  ','N1  ','N6  ','N7  ',5,'C4  ','C2  ',
     1  'N1  ','C4  ','C','w:w','w:w',0)

c       ADE-ADE (whC) HB betn. N7A-N6A & C8A-N1A
        call findpair(i,1,'C2  ','N7  ','N6  ','N7  ',1,'C4  ','C8  ',
     1  'N1  ','C4  ','C','h:w','w:h',0)

c       ADE-N6G (whC) HB betn. N7A-N6A & C8A-N1A
        call findpair(i,1,'C2  ','N7  ','N6  ','N7  ',5,'C4  ','C8  ',
     1  'N1  ','C4  ','C','h:w','w:h',0)

c       N6G-ADE (whC) HB betn. N7A-N6A & C8A-N1A
        call findpair(i,5,'C2  ','N7  ','N6  ','N7  ',1,'C4  ','C8  ',
     1  'N1  ','C4  ','C','h:w','w:h',0)

c       N6G-N6G (whC) HB betn. N7A-N6A & C8A-N1A
        call findpair(i,5,'C2  ','N7  ','N6  ','N7  ',5,'C4  ','C8  ',
     1  'N1  ','C4  ','C','h:w','w:h',0)

c       ADE-ADE (wsC) HB betn. N1A-C2A & C2A-N3A ! Removed as both CHO
c        call findpair(i,1,'C5  ','N1  ','C2  ','C6  ',1,'C4  ','C2  ',
c     1  'N3  ','C5  ','C','w:s','s:w',0)        ! bhatta May 2019

c       ADE-ADE (wsT) HB betn. N6A-N3A & N1A-C2A 
        call findpair(i,1,'N7  ','N6  ','N3  ','C5  ',1,'C4  ','N1  ',
     1  'C2  ','C6  ','T','w:s','s:w',0)

c       N6G-ADE (wsT) HB betn. N6A-N3A & N1A-C2A 
        call findpair(i,5,'N7  ','N6  ','N3  ','C5  ',1,'C4  ','N1  ',
     1  'C2  ','C6  ','T','w:s','s:w',0)

c       ADE-ADE (hhC) HB betn. N7A-N6A & C8A-N7A
        call findpair(i,1,'C4  ','N7  ','N6  ','N1  ',1,'C4  ','C8  ',
     1  'N7  ','N9  ','C','h:h','h:h',0)

c       ADE-N6G (hhC) HB betn. N7A-N6A & C8A-N7A
        call findpair(i,1,'C4  ','N7  ','N6  ','N1  ',5,'C4  ','C8  ',
     1  'N7  ','N9  ','C','h:h','h:h',0)

c       N6G-N6G (hhC) HB betn. N7A-N6A & C8A-N7A
        call findpair(i,5,'C4  ','N7  ','N6  ','N1  ',5,'C4  ','C8  ',
     1  'N7  ','N9  ','C','h:h','h:h',0)

c       ADE-ADE (hsT) HB betn. N6A-N3A & N7A-C2A
        call findpair(i,1,'N1  ','N6  ','N3  ','C6  ',1,'N9  ','N7  ',
     1  'C2  ','C5  ','T','h:s','s:h',0)

c       N6G-ADE (hsT) HB betn. N6A-N3A & N7A-C2A
        call findpair(i,5,'N1  ','N6  ','N3  ','C6  ',1,'N9  ','N7  ',
     1  'C2  ','C5  ','T','h:s','s:h',0)

c       ADE-ADE (ssT) HB betn. N3A-C2A & C2A-N3A ! Removed as both CHO
c        call findpair(i,1,'C5  ','N3  ','C2  ','C6  ',1,'C6  ','C2  ',
c     1  'N3  ','C5  ','T','s:s','s:s',0)       ! bhatta May 2019

c*******************************ADE-GUA********************************

c       ADE-GUA (wwT) HB betn. N1A-N1G & C2A-O6G
        call findpair(i,1,'C5  ','N1  ','N1  ','C4  ',2,'C4  ','C2  ',
     1  'O6  ','N7  ','T','w:w','w:w',0)

c       ADE-QUO (wwT) HB betn. N1A-N1G & C2A-O6G
        call findpair(i,1,'C5  ','N1  ','N1  ','C4  ',7,'C4  ','C2  ',
     1  'O6  ','C7  ','T','w:w','w:w',0)

c       ADE-GUA (wsC) HB betn. N1A-N2G & C2A-N3G & N3A-C1*G Involves wsC & ssC
c        call findhb3(i,1,2,'C5  ','N1  ','N2  ','N1  ','C4  ','C2  ',
c     1  'N3  ','C5  ','C5  ','N3  ','C1* ','C8  ','C','w:s','s:w')

c       ADE-GUA (wsC) HB betn. N1A-N2G & C2A-N3G
        call findpair(i,1,'C5  ','N1  ','N2  ','N1  ',2,'C4  ','C2  ',
     1  'N3  ','C5  ','C','w:s','s:w',0)

c       ADE-QUO (wsC) HB betn. N1A-N2G & C2A-N3G
        call findpair(i,1,'C5  ','N1  ','N2  ','N1  ',7,'C4  ','C2  ',
     1  'N3  ','C5  ','C','w:s','s:w',0)

c       ADE-N6G (wsC) HB betn. N1A-N2G & C2A-N3G
        call findpair(i,1,'C5  ','N1  ','N2  ','N1  ',5,'C4  ','C2  ',
     1  'N3  ','C5  ','C','w:s','s:w',0)

c       ADE-GUA (hwT) HB betn. N7A-N1G & C8A-O6G
        call findpair(i,1,'N3  ','N7  ','N1  ','C4  ',2,'C4  ','C8  ',
     1  'O6  ','N7  ','T','h:w','w:h',0)

c       ADE-QUO (hwT) HB betn. N7A-N1G & C8A-O6G
        call findpair(i,1,'N3  ','N7  ','N1  ','C4  ',7,'C4  ','C8  ',
     1  'O6  ','C7  ','T','h:w','w:h',0)

c       N6G-GUA (hwT) HB betn. N7A-N1G & C8A-O6G
        call findpair(i,5,'N3  ','N7  ','N1  ','C4  ',2,'C4  ','C8  ',
     1  'O6  ','N7  ','T','h:w','w:h',0)

c       N6G-QUO (hwT) HB betn. N7A-N1G & C8A-O6G
        call findpair(i,5,'N3  ','N7  ','N1  ','C4  ',7,'C4  ','C8  ',
     1  'O6  ','C7  ','T','h:w','w:h',0)

c       ADE-GUA (hhC) HB betn. N6A-N7G & N7A-C8G
        call findpair(i,1,'N1  ','N6  ','N7  ','C4  ',2,'C4  ','N7  ',
     1  'C8  ','C4  ','C','h:h','h:h',0)

c       ADE-GUA (hsC) HB betn. N7A-N2G & C8A-N3G
        call findpair(i,1,'C5  ','N7  ','N2  ','N1  ',2,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:s','s:h',0)

c       ADE-QUO (hsC) HB betn. N7A-N2G & C8A-N3G
        call findpair(i,1,'C5  ','N7  ','N2  ','N1  ',7,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:s','s:h',0)

c       N6G-GUA (hsC) HB betn. N7A-N2G & C8A-N3G
        call findpair(i,5,'C5  ','N7  ','N2  ','N1  ',2,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:s','s:h',0)

c       N6G-QUO (hsC) HB betn. N7A-N2G & C8A-N3G
        call findpair(i,5,'C5  ','N7  ','N2  ','N1  ',7,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:s','s:h',0)

c       ADE-N6G (hsC) HB betn. N7A-N2G & C8A-N3G
        call findpair(i,1,'C5  ','N7  ','N2  ','N1  ',5,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:s','s:h',0)

c       N6G-N6G (hsC) HB betn. N7A-N2G & C8A-N3G
        call findpair(i,5,'C5  ','N7  ','N2  ','N1  ',5,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:s','s:h',0)

c       ADE-GUA (swC) HB betn. C2A-O6G & N3A-N1G
        call findpair(i,1,'C6  ','C2  ','O6  ','N7  ',2,'C5  ','N3  ',
     1  'N1  ','C4  ','C','s:w','w:s',0) !Corrected N3 to N7,bhatta 2019

c       ADE-GUA (swT) HB betn. N3A-N1G & C2A-O6G
        call findpair(i,1,'C5  ','N3  ','N1  ','N3  ',2,'C6  ','C2  ',
     1  'O6  ','N7  ','C','s:w','w:s',0)

c       ADE-QUO (swT) HB betn. N3A-N1G & C2A-O6G
        call findpair(i,1,'C5  ','N3  ','N1  ','N3  ',7,'C6  ','C2  ',
     1  'O6  ','C7  ','C','s:w','w:s',0)

c       ADE-GUA (ssT) HB betn. N3A-N2G & C2A-N3G 
        call findpair(i,1,'C5  ','N3  ','N2  ','N1  ',2,'C5  ','C2  ',
     1  'N3  ','C5  ','T','s:s','s:s',0)

c       ADE-QUO (ssT) HB betn. N3A-N2G & C2A-N3G
        call findpair(i,1,'C5  ','N3  ','N2  ','N1  ',7,'C5  ','C2  ',
     1  'N3  ','C5  ','T','s:s','s:s',0)

c       ADE-N6G (ssT) HB betn. N3A-N2G & C2A-N3G
        call findpair(i,1,'C5  ','N3  ','N2  ','N1  ',5,'C5  ','C2  ',
     1  'N3  ','C5  ','T','s:s','s:s',0)

c       ADE-GUA (zhC) HB betn. N3A-N7G & C2A-O6G
        call findpair(i,1,'C6  ','N3  ','N7  ','N9  ',2,'C6  ','C2  ',
     1  'O6  ','N1  ','C','z:h','h:z',0)

c       ADE-GUA (zhT) HB betn. N3A-O6G & C2A-N7G
        call findpair(i,1,'C5  ','N3  ','O6  ','C1  ',2,'C6  ','C2  ',
     1  'N7  ','N9  ','T','z:h','h:z',0)

c*******************************ADE-CYT********************************

c       ADE-CYT (wwC) HB betn. N1A-N4C & C2A-N3C
        call findpair(i,1,'C5  ','N1  ','N4  ','C5  ',3,'C8  ','C2  ',
     1  'N3  ','C6  ','C','w:w','w:w',0)

c       ADE-CYT (hwC) HB betn. N7A-N4C & C8A-N3C
        call findpair(i,1,'C5  ','N7  ','N4  ','C5  ',3,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:w','w:h',0)

c       N6G-CYT (hwC) HB betn. N7A-N4C & C8A-N3C
        call findpair(i,5,'C5  ','N7  ','N4  ','C5  ',3,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:w','w:h',0)

c       ADE-CYT (swT) HB betn. N3A-N4C & C2A-N3C
        call findpair(i,1,'C5  ','N3  ','N4  ','N3  ',3,'C6  ','C2  ',
     1  'N3  ','C6  ','T','s:w','w:s',0)

c*******************************ADE-URA********************************

c       ADE-URA (whC) HB betn. N6A-O4U & N1A-C5U
        call findpair(i,1,'N7  ','N6  ','O4  ','N3  ',4,'C8  ','N1  ',
     1  'C5  ','C2  ','C','w:h','h:w',0)

c       N6G-URA (whC) HB betn. N6A-O4U & N1A-C5U
        call findpair(i,5,'N7  ','N6  ','O4  ','N3  ',4,'C8  ','N1  ',
     1  'C5  ','C2  ','C','w:h','h:w',0)

c       ADE-URA (whT) HB betn. N1A-C5U & C2A-O4U
        call findpair(i,1,'C8  ','N1  ','C2  ','C5  ',4,'C4  ','C2  ',
     1  'O4  ','N3  ','T','w:h','h:w',0)

c       ADE-PSU (whT) HB betn. N1A-N1U & C2A-O2U
        call findpair(i,1,'C5  ','N1  ','N1  ','C4  ',6,'C4  ','C2  ',
     1  'O2  ','N3  ','T','w:h','h:w',0)

c       ADE-URA (hhC) HB betn. N6A-O4U & N7A-C5U
        call findpair(i,1,'N1  ','N6  ','O4  ','N3  ',4,'C4  ','N7  ',
     1  'C5  ','N1  ','C','h:h','h:h',0)

c       N6G-URA (hhC) HB betn. N6A-O4U & N7A-C5U
        call findpair(i,5,'N1  ','N6  ','O4  ','N3  ',4,'C4  ','N7  ',
     1  'C5  ','N1  ','C','h:h','h:h',0)

c       ADE-URA (hhT) HB betn. N7A-C5U & C8A-O4U
        call findpair(i,1,'C4  ','N7  ','C5  ','N1  ',4,'C4  ','C8  ',
     1  'O4  ','N3  ','T','h:h','h:h',0)

c       N6G-URA (hhT) HB betn. N7A-C5U & C8A-O4U
        call findpair(i,5,'C4  ','N7  ','C5  ','N1  ',4,'C4  ','C8  ',
     1  'O4  ','N3  ','T','h:h','h:h',0)

c       ADE-URA (swC) HB betn. N3A-N3U & C2A-O4U
        call findpair(i,1,'C5  ','N3  ','N3  ','C6  ',4,'C6  ','C2  ',
     1  'O4  ','C5  ','C','s:w','w:s',0) ! Modified C2U to C6U, bhatta 2019

c       ADE-PSU (swC) HB betn. N3A-N3U & C2A-O2U
        call findpair(i,1,'C5  ','N3  ','N3  ','C5  ',6,'C6  ','C2  ',
     1  'O2  ','N1  ','C','s:w','w:s',0)

c       ADE-URA (swT) HB betn. N3A-N3U & C2A-O2U
        call findpair(i,1,'C5  ','N3  ','N3  ','C5  ',4,'C6  ','C2  ',
     1  'O2  ','C1* ','T','s:w','w:s',0)

c       ADE-PSU (swT) HB betn. N3A-N3U & C2A-O4U
        call findpair(i,1,'C5  ','N3  ','N3  ','C2  ',6,'C6  ','C2  ',
     1  'O4  ','C1* ','T','s:w','w:s',0)

c       ADE-URA (shC) HB betn. N3A-C5U & C2A-O4U
        call findpair(i,1,'C5  ','N3  ','C5  ','N1  ',4,'C6  ','C2  ',
     1  'O4  ','N3  ','C','s:h','h:s',0)

c       ADE-PSU (shC) HB betn. N3A-N1U & C2A-O2U
        call findpair(i,1,'C5  ','N3  ','N1  ','C5  ',6,'C6  ','C2  ',
     1  'O2  ','N3  ','C','s:h','h:s',0)

c*******************************GUA-GUA********************************

c       GUA-GUA (hhT) HB betn. N7G-C8G & C8G-N7G 
        call findpair(i,2,'C2  ','N7  ','C8  ','C4  ',2,'C4  ','C8  ',
     1  'N7  ','C2  ','T','h:h','h:h',0)

c       GUA-GUA (hsC) HB betn. N7G-N2G & C8G-N3G
        call findpair(i,2,'C4  ','N7  ','N2  ','N1  ',2,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:s','s:h',0)

c       GUA-QUO (hsC) HB betn. N7G-N2G & C8G-N3G
        call findpair(i,2,'C4  ','N7  ','N2  ','N1  ',7,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:s','s:h',0)

c       GUA-N6G (hsC) HB betn. N7G-N2G & C8G-N3G
        call findpair(i,2,'C4  ','N7  ','N2  ','N1  ',5,'C4  ','C8  ',
     1  'N3  ','C6  ','C','h:s','s:h',0)

c*******************************GUA-CYT********************************

c       GUA-CYT (hhC) HB betn. O6G-N4C & N7G-C5C
        call findpair(i,2,'N1  ','O6  ','N4  ','N3  ',3,'N9  ','N7  ',
     1  'C5  ','C2  ','C','h:h','h:h',0)

c       GUA-CYT (hhT) HB betn. N7G-N4C & O6G-C5C 
        call findpair(i,2,'N9  ','N7  ','N4  ','N3  ',3,'N1  ','O6  ',
     1  'C5  ','N1  ','T','h:h','h:h',0)

c*******************************GUA-URA********************************

c       GUA-URA (whT) HB betn. N1G-O4U & O6G-C5U
        call findpair(i,2,'C4  ','N1  ','O4  ','N3  ',4,'N7  ','O6  ',
     1  'C5  ','N1  ','T','w:h','h:w',0)

c       QUO-URA (whT) HB betn. N1G-O4U & O6G-C5U
        call findpair(i,7,'C4  ','N1  ','O4  ','N3  ',4,'C7  ','O6  ',
     1  'C5  ','N1  ','T','w:h','h:w',0)

c       GUA-URA (hwC) HB betn. N7G-N3U & C8G-O2U
        call findpair(i,2,'C4  ','N7  ','N3  ','C6  ',4,'N9  ','C8  ',
     1  'O2  ','N1  ','C','h:w','w:h',0)

c       GUA-PSU (hwC) HB betn. N7G-N3U & C8G-O4U
        call findpair(i,2,'C5  ','N7  ','N3  ','C6  ',6,'C4  ','C8  ',
     1  'O4  ','C1* ','C','h:w','w:h',0)

c       GUA-URA (hwT) HB betn. N7G-N3U & C8G-O4U
        call findpair(i,2,'C5  ','N7  ','N3  ','C6  ',4,'C4  ','C8  ',
     1  'O4  ','C5  ','T','h:w','w:h',0)

c       GUA-PSU (hwT) HB betn. N7G-N3U & C8G-O2U
        call findpair(i,2,'C4  ','N7  ','N3  ','C6  ',6,'N9  ','C8  ',
     1  'O2  ','N1  ','T','h:w','w:h',0)

c       GUA-URA (hhT) HB betn. N7G-C5U & C8G-O4U
        call findpair(i,2,'C2  ','N7  ','C5  ','N1  ',4,'C4  ','C8  ',
     1  'O4  ','N3  ','T','h:h','h:h',0)

c       GUA-URA (shC) HB betn. N2G-O4U & N3G-C5U 
        call findpair(i,2,'N1  ','N2  ','O4  ','N3  ',4,'C6  ','N3  ',
     1  'C5  ','C2  ','C','s:h','h:s',0)

c       QUO-URA (shC) HB betn. N2G-O4U & N3G-C5U 
        call findpair(i,7,'N1  ','N2  ','O4  ','N3  ',4,'C6  ','N3  ',
     1  'C5  ','C2  ','C','s:h','h:s',0)

c       N6G-URA (shC) HB betn. N2G-O4U & N3G-C5U 
        call findpair(i,5,'N1  ','N2  ','O4  ','N3  ',4,'C6  ','N3  ',
     1  'C5  ','C2  ','C','s:h','h:s',0)

c*******************************CYT-CYT********************************

c       CYT-CYT (whC) HB betn. N3C-N4C & O2C-C5C 
        call findpair(i,3,'C6  ','N3  ','N4  ','N3  ',3,'N1  ','O2  ',
     1  'C5  ','C2  ','C','w:h','h:w',0)

c       CYT-CYT (whT) HB betn. O2C-N4C & N3C-C5C
        call findpair(i,3,'N1  ','O2  ','N4  ','N3  ',3,'C6  ','N3  ',
     1  'C5  ','C2  ','T','w:h','h:w',0)

c*******************************CYT-URA********************************

c       CYT-URA (whC) HB betn. N4C-O4U & N3C-C5U 
        call findpair(i,3,'C5  ','N4  ','O4  ','N3  ',4,'C6  ','N3  ',
     1  'C5  ','C2  ','C','w:h','h:w',0)

c*******************************URA-URA********************************

c       URA-URA (whC) HB betn. N3U-O4U & O2U-C5U 
        call findpair(i,4,'C5  ','N3  ','O4  ','N3  ',4,'N1  ','O2  ',
     1  'C5  ','C2  ','C','w:h','h:w',0)

c       PSU-URA (whC) HB betn. N3U-O4U & O4U-C5U
        call findpair(i,6,'N1  ','N3  ','O4  ','N3  ',4,'C5  ','O4  ',
     1  'C5  ','N1  ','C','w:h','h:w',0)

c       URA-URA (whT) HB betn. N3U-O4U & O4U-C5U
        call findpair(i,4,'N1  ','N3  ','O4  ','N3  ',4,'C5  ','O4  ',
     1  'C5  ','N1  ','T','w:h','h:w',0)

c       PSU-URA (whT) HB betn. N3U-O4U & O2U-C5U 
        call findpair(i,6,'C5  ','N3  ','O4  ','N3  ',4,'N1  ','O2  ',
     1  'C5  ','C2  ','T','w:h','h:w',0)

c       URA-URA (hhT) HB betn. O4U-C5U & C5U-O4U
        call findpair(i,4,'N3  ','O4  ','C5  ','N1  ',4,'N1  ','C5  ',
     1  'O4  ','N3  ','T','h:h','h:h',0)

c       URA-PSU (hhT) HB betn. O4U-N1U & C5U-O2U
        call findpair(i,4,'C2  ','O4  ','N1  ','C5  ',6,'N1  ','C5  ',
     1  'O2  ','N3  ','T','h:h','h:h',0)

c**********************************************************************

	endif

c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%
c	---  Base pairs having C-H...O H-bonds involving sugar O2* ---
c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if(chohb.eq.1.and.sugarbp.eq.1) then

c*******************************ADE-ADE********************************

c       ADE-ADE (ssC) HB betn. N3A-C2A & C1*A-N3A
        call findpair(i,1,'C6  ','N3  ','C2  ','C5  ',1,'C8  ','C1* ',
     1  'N3  ','N7  ','C','s:s','s:s',1)

c*******************************ADE-GUA********************************

c       ADE-GUA (wsT) HB betn. C2A-N3G & N1A-C1*G
c       call findpair(i,1,'C5  ','C2  ','N3  ','C5  ',2,'C5  ','N1  ',
c    1  'C1* ','C8  ','T','s:s','s:s',1)    Sukanya
        call findpair(i,2,'C5  ','N3  ','C2  ','C5  ',1,'C8  ','C1* ',
     1  'N1  ','C5  ','T','s:w','w:s',1)

c       ADE-GUA (ssC) HB betn. C2A-N3G & N3A-C1*G
c       call findpair(i,1,'C4  ','C2  ','N3  ','C5  ',2,'C5  ','N3  ',
c    1  'C1* ','C8  ','C','s:s','s:s',1)    Sukanya
        call findpair(i,2,'C5  ','N3  ','C2  ','C4  ',1,'C8  ','C1* ',
     1  'N3  ','C5  ','C','s:s','s:s',1)

c*******************************ADE-CYT********************************

c       ADE-CYT (wsT) HB betn. C2A-O2C & N6A-C1*C or N1A-C1*C
        call findhb3(i,3,1,'N3  ','O2  ','C2  ','C4  ','C6  ','C1* ',
     1  'N6  ','N7  ','C6  ','C1* ','N1  ','C5  ','T','s:w','w:s')

c       ADE-CYT (hsT) HB betn. C8A-O2C & N7A-C1*C
c       call findpair(i,1,'C4  ','C8  ','O2  ','C4  ',3,'C4  ','N7  ',
c    1  'C1* ','C6  ','T','h:s','s:h',1)    Sukanya
        call findpair(i,3,'C4  ','O2  ','C8  ','C4  ',1,'C6  ','C1* ',
     1  'N7  ','C4  ','T','s:h','h:s',1)

c       N6G-CYT (hsT) HB betn. C8A-O2C & N7A-C1*C
c       call findpair(i,5,'C4  ','C8  ','O2  ','C4  ',3,'C4  ','N7  ',
c    1  'C1* ','C6  ','T','h:s','s:h',1)    Sukanya
        call findpair(i,3,'C4  ','O2  ','C8  ','C4  ',5,'C6  ','C1* ',
     1  'N7  ','C4  ','T','s:h','h:s',1)

c       ADE-CYT (shC) HB betn. N3A-N4C & C1*A-C5C
        call findpair(i,1,'C5  ','N3  ','N4  ','N3  ',3,'C8  ','C1* ',
     1  'C5  ','C2  ','C','s:h','h:s',1)

c       ADE-CYT (shT) HB betn. N3A-C5C & C1*A-N4C
        call findpair(i,1,'C5  ','N3  ','C5  ','C2  ',3,'C8  ','C1* ',
     1  'N3  ','N4  ','T','s:h','h:s',1)

c       ADE-CYT (ssC) HB betn. C2A-O2C & N3A-C1*C
c       call findpair(i,1,'C5  ','C2  ','O2  ','C5  ',3,'C4  ','N3  ',
c    1  'C1* ','C6  ','C','s:s','s:s',1)    Sukanya
        call findpair(i,3,'C5  ','O2  ','C2  ','C5  ',1,'C6  ','C1* ',
     1  'N3  ','C4  ','C','s:s','s:s',1)

c*******************************ADE-URA********************************

c       ADE-PSU (wsC) HB betn. N6A-O4U & C2A-C1*U
c       call findpair(i,1,'N7  ','N6  ','O4  ','N1  ',6,'N9  ','C2  ',
c    1  'C1* ','C6  ','C','w:s','s:w',1)    Sukanya
        call findpair(i,6,'N1  ','O4  ','N6  ','N7  ',1,'C6  ','C1* ',
     1  'C2  ','N9  ','C','s:w','w:s',1)

c       ADE-URA (wsT) HB betn. C2A-O2U & N6A-C1*U
        call findpair(i,4,'N3  ','O2  ','C2  ','C5  ',1,'C6  ','C1* ',
     1  'N6  ','N7  ','T','s:w','w:s',1)

c       ADE-URA (wsT) HB betn. C2A-O2U & N1A-C1*U
        call findpair(i,4,'N3  ','O2  ','C2  ','C5  ',1,'C6  ','C1* ',
     1  'N1  ','C5  ','T','s:w','w:s',1)

c       ADE-PSU (wsT) HB betn. C2A-O4U & N6A-C1*U
c       call findpair(i,1,'N9  ','C2  ','O4  ','N1  ',6,'N7  ','N6  ',
c    1  'C1* ','C6  ','T','w:s','s:w',1)    Sukanya
        call findpair(i,6,'N1  ','O4  ','C2  ','N9  ',1,'C6  ','C1* ',
     1  'N6  ','N7  ','T','s:w','w:s',1)

c       ADE-URA (hsT) HB betn. C8A-O2U & N6A-C1*U
        call findpair(i,4,'C5  ','O2  ','C8  ','C4  ',1,'C6  ','C1* ',
     1  'N6  ','N1  ','T','s:h','h:s',1)

c       N6G-URA (hsC) HB betn. N6A-O2U & C8A-C1*U
        call findpair(i,4,'N3  ','O2  ','N6  ','N1  ',5,'C6  ','C1* ',
     1  'C8  ','C1* ','C','s:h','h:s',1)

c       ADE-URA (shT) HB betn. N3A-C5U & C1*A-O4U
        call findpair(i,1,'C5  ','N3  ','C5  ','C2  ',4,'C8  ','C1* ',
     1  'O4  ','N3  ','T','s:h','h:s',1)

c       ADE-URA (ssC) HB betn. C2A-O2U & N3A-C1*U
c       call findpair(i,1,'C5  ','C2  ','O2  ','C5  ',4,'N7  ','N3  ',
c    1  'C1* ','C6  ','C','s:s','s:s',1)    Sukanya
        call findpair(i,4,'C5  ','O2  ','C2  ','C5  ',1,'C6  ','C1* ',
     1  'N3  ','N7  ','C','s:s','s:s',1)

c*******************************GUA-GUA********************************

c       GUA-GUA (hsT) HB betn. N3G-C8G & C1*G-N7G
c        call findpair(i,2,'C5  ','N3  ','C8  ','C4  ',2,'C8  ','C1* ',
c     1  'N7  ','N3  ','T','h:s','s:h',1)  Removed.  Structures obtained are 
c  too bad with bad NUPARM parameters.  bhatta March 4, 2010.

c       GUA-QUO (hsT) HB betn. N3G-C8G & C1*G-N7G
        call findpair(i,7,'C5  ','N3  ','C8  ','C4  ',2,'C8  ','C1* ',
     1  'N7  ','N3  ','T','s:h','h:s',1)

c       N6G-GUA (hsT) HB betn. N3G-C8G & C1*G-N7G
        call findpair(i,5,'C5  ','N3  ','C8  ','C4  ',2,'C8  ','C1* ',
     1  'N7  ','N3  ','T','h:s','s:h',1)

c*******************************GUA-CYT********************************

c       GUA-CYT (hsT) HB betn. C8G-O2C & N7G-C1*C
c       call findpair(i,2,'C4  ','C8  ','O2  ','C5  ',3,'N1  ','N7  ',
c    1  'C1* ','C6  ','T','h:s','s:h',1)    Sukanya
        call findpair(i,3,'C5  ','O2  ','C8  ','C4  ',2,'C6  ','C1* ',
     1  'N7  ','N1  ','T','s:h','h:s',1)

c*******************************GUA-URA********************************

c       GUA-URA (hsT) HB betn. C8G-O2U & N7G-C1*U
c       call findpair(i,2,'C4  ','C8  ','O2  ','N3  ',4,'C2  ','N7  ',
c    1  'C1* ','C6  ','T','h:s','s:h',1)    Sukanya
        call findpair(i,4,'N3  ','O2  ','C8  ','C4  ',4,'C6  ','C1* ',
     1  'N7  ','C2  ','T','s:h','h:s',1)

c*******************************CYT-CYT********************************

c       CYT-CYT (shC) HB betn. O2C-N4C & C1*C-C5C
        call findpair(i,3,'N3  ','O2  ','N4  ','N3  ',3,'C6  ','C1* ',
     1  'C5  ','C1* ','C','s:h','h:s',1)

c       CYT-CYT (shT) HB betn. O2C-C5C & C1*C-N4C
        call findpair(i,3,'N3  ','O2  ','C5  ','N1  ',3,'C6  ','C1* ',
     1  'N4  ','N3  ','T','s:h','h:s',1)

c*******************************CYT-URA********************************

c       CYT-URA (hsC) HB betn. N4C-O2U & C5C-C1*U
        call findpair(i,4,'N3  ','O2  ','N4  ','N3  ',3,'C6  ','C1* ',
     1  'C5  ','C1* ','C','s:h','h:s',1)

c       CYT-URA (hsT) HB betn. C5C-O2U & N4C-C1*U
        call findpair(i,4,'N3  ','O2  ','C5  ','N1  ',3,'C6  ','C1* ',
     1  'N4  ','N3  ','T','s:h','h:s',1)

c       CYT-PSU (hsT) HB betn. C5C-O4U & N4C-C1*U
c       call findpair(i,3,'C2  ','C5  ','O4  ','N1  ',6,'N3  ','N4  ',
c    1  'C1* ','C6  ','T','h:s','s:h',1)    Sukanya
        call findpair(i,6,'N1  ','O4  ','C5  ','C2  ',3,'C6  ','C1* ',
     1  'N4  ','N3  ','T','s:h','h:s',1)

c       CYT-URA (shT) HB betn. O2C-C5U & C1*C-O4U
        call findpair(i,3,'C5  ','O2  ','C5  ','C2  ',4,'C6  ','C1* ',
     1  'O4  ','N3  ','T','s:h','h:s',1)

c*******************************URA-URA********************************

c       URA-URA (hsT) HB betn. C5U-O2U & O4U-C1*U
c       call findpair(i,4,'C2  ','C5  ','O2  ','C5  ',4,'N3  ','O4  ',
c    1  'C1* ','C6  ','T','h:s','s:h',1)    Sukanya
        call findpair(i,4,'C5  ','O2  ','C5  ','C2  ',4,'C6  ','C1* ',
     1  'O4  ','N3  ','T','s:h','h:s',1)

c       URA-PSU (hsT) HB betn. C5U-O4U & O4U-C1*U
c       call findpair(i,4,'C2  ','C5  ','O4  ','C2  ',6,'N3  ','O4  ',
c    1  'C1* ','C6  ','T','h:s','s:h',1)    Sukanya
        call findpair(i,6,'C2  ','O4  ','C5  ','C2  ',4,'C6  ','C1* ',
     1  'N3  ','O4  ','T','s:h','h:s',1)

c**********************************************************************

        endif

c *******************************************************************
c        Base pairs stabilized by N-H...O/N H-bonds through sugar O2*
c *******************************************************************

	if(sugarbp.eq.1) then

c*******************************ADE-ADE********************************

c       ADE-ADE (HSC) HB betn. N6A-N3A & N7A-C1*A
c       call findpair(i,1,'N1  ','N6  ','N3  ','C5  ',1,'C4  ','N7  ',
c    1  'C1* ','C8  ','C','H:S','S:H',1)    Sukanya
        call findpair(i,1,'C5  ','N3  ','N6  ','N1  ',1,'C8  ','C1* ',
     1  'N7  ','C4  ','C','S:H','H:S',1)

c       ADE-ADE (WSC) HB betn. N6A-N3A & N1A-C1*A
c       call findpair(i,1,'N7  ','N6  ','N3  ','C5  ',1,'C4  ','N1  ',
c    1  'C1* ','C8  ','C','W:S','S:W',1)     ! Considering SAM in 3IQN.pdb Sukanya
        call findpair(i,1,'C5  ','N3  ','N6  ','N7  ',1,'C8  ','C1* ',
     1  'N1  ','C4  ','C','S:W','W:S',1)     ! Considering SAM in 3IQN.pdb

c       N6G-ADE (WSC) HB betn. N6A-N3A & N1A-C1*A
c       call findpair(i,5,'N7  ','N6  ','N3  ','C5  ',1,'C4  ','N1  ',
c    1  'C1* ','C8  ','C','W:S','S:W',1)    Sukanya 
        call findpair(i,5,'C5  ','N3  ','N6  ','N7  ',1,'C8  ','C1* ',
     1  'N1  ','C4  ','C','S:W','W:S',1)

c       N6G-ADE (HSC) HB betn. N6A-N3A & N7A-C1*A
c       call findpair(i,5,'N1  ','N6  ','N3  ','C5  ',1,'C4  ','N7  ',
c    1  'C1* ','C8  ','C','H:S','S:H',1)    Sukanya
        call findpair(i,1,'C5  ','N3  ','N6  ','N1  ',5,'C8  ','C1* ',
     1  'N7  ','C4  ','C','S:H','H:S',1)

c*******************************ADE-GUA********************************

c       ADE-GUA (SWC) HB betn. N3A-N1G & C1*A-N2G
        call findpair(i,1,'C5  ','N3  ','N1  ','C4  ',2,'C8  ','C1* ',
     1  'N2  ','N3  ','C','S:W','W:S',1)

c       ADE-GUA (SSC) HB betn. N3A-N2G & C1*A-N3G
        call findpair(i,1,'C5  ','N3  ','N2  ','N1  ',2,'C8  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c*******************************ADE-CYT********************************

c       ADE-CYT (WSC) HB betn. N6A-O2C & N1A-C1*C
        call findpair(i,3,'N3  ','O2  ','N6  ','N7  ',1,'C6  ','C1* ',
     1  'N1  ','C4  ','C','S:W','W:S',1)

c       ADE-CYT (HSC) HB betn. N6A-O2C & N7A-C1*C
        call findpair(i,3,'N3  ','O2  ','N6  ','N1  ',1,'C6  ','C1* ',
     1  'N7  ','C4  ','C','S:H','H:S',1)

c       N6G-CYT (HSC) HB betn. N6A-O2C & N7A-C1*C
        call findpair(i,3,'N3  ','O2  ','N6  ','N1  ',5,'C6  ','C1* ',
     1  'N7  ','C4  ','C','S:H','H:S',1)

c       ADE-CYT (SWC) HB betn. N3A-N4C & C1*A-N3C
        call findpair(i,1,'C5  ','N3  ','N4  ','C5  ',3,'C8  ','C1* ',
     1  'N3  ','C6  ','C','S:W','W:S',1)

c       ADE-CYT (SST) HB betn. N3A-C1*C & C1*A-O2C Removed bhatta June 2019
C        call findpair(i,1,'C5  ','N3  ','C1* ','C6  ',3,'C8  ','C1* ',
C     1  'O2  ','C5  ','T','S:S','S:S',1)

c*******************************ADE-URA********************************

c       ADE-URA (WSC) HB betn. N6A-O2U & N1A-C1*U
        call findpair(i,4,'C5  ','O2  ','N6  ','N7  ',1,'C6  ','C1* ',
     1  'N1  ','C4  ','C','S:W','W:S',1)

c       ADE-URA (HSC) HB betn. N6A-O2U & N7A-C1*U
c       call findpair(i,1,'N1  ','N6  ','O2  ','C5  ',4,'C4  ','N7  ',
c    1  'C1* ','C6  ','C','H:S','S:H',1)    Sukanya
        call findpair(i,4,'C5  ','O2  ','N6  ','N1  ',1,'C6  ','C1* ',
     1  'N7  ','C4  ','C','S:H','H:S',1)

c       N6G-URA (HSC) HB betn. N6A-O2U & N7A-C1*U
c       call findpair(i,5,'N1  ','N6  ','O2  ','C5  ',4,'C4  ','N7  ',
c    1  'C1* ','C6  ','C','H:S','S:H',1)    Sukanya
        call findpair(i,4,'C5  ','O2  ','N6  ','N1  ',5,'C6  ','C1* ',
     1  'N7  ','C4  ','C','S:H','H:S',1)

c       ADE-PSU (HSC) HB betn. N6A-O4U & N7A-C1*U
c       call findpair(i,1,'N1  ','N6  ','O4  ','N1  ',6,'C4  ','N7  ',
c    1  'C1* ','C6  ','C','H:S','S:H',1)    Sukanya
        call findpair(i,6,'N1  ','O4  ','N6  ','N1  ',1,'C6  ','C1* ',
     1  'N7  ','C4  ','C','S:H','H:S',1)

c       N6G-PSU (HSC) HB betn. N6A-O4U & N7A-C1*U
c       call findpair(i,5,'N1  ','N6  ','O4  ','N1  ',6,'C4  ','N7  ',
c    1  'C1* ','C6  ','C','H:S','S:H',1)    Sukanya
        call findpair(i,6,'N1  ','O4  ','N6  ','N1  ',5,'C6  ','C1* ',
     1  'N7  ','C4  ','C','S:H','H:S',1)

c       ADE-PSU (SHT) HB betn. N3A-N1U & C1*A-O2U
        call findpair(i,1,'C5  ','N3  ','N1  ','C4  ',6,'C8  ','C1* ',
     1  'O2  ','N3  ','T','S:H','H:S',1)

c       ADE-URA (SST) HB betn. N3A-C1*U & C1*A-O2U
        call findpair(i,1,'C5  ','N3  ','C1* ','C6  ',4,'C8  ','C1* ',
     1  'O2  ','C5  ','T','S:S','S:S',1)

c*******************************GUA-GUA********************************

c       GUA-GUA (WST) HB betn. N1G-N3G & O6G-C1*G
c       call findpair(i,2,'C4  ','N1  ','N3  ','C6  ',2,'N7  ','O6  ',
c    1  'C1* ','C8  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,2,'C6  ','N3  ','N1  ','C4  ',2,'C8  ','C1* ',
     1  'O6  ','N7  ','T','S:W','W:S',1)

c       QUO-GUA (WST) HB betn. N1G-N3G & O6G-C1*G
c       call findpair(i,7,'C4  ','N1  ','N3  ','C6  ',2,'C7  ','O6  ',
c    1  'C1* ','C8  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,2,'C6  ','N3  ','N1  ','C4  ',7,'C8  ','C1* ',
     1  'O6  ','C7  ','T','S:W','W:S',1)

c       GUA-QUO (WST) HB betn. N1G-N3G & O6G-C1*G
c       call findpair(i,2,'C4  ','N1  ','N3  ','C6  ',7,'N7  ','O6  ',
c    1  'C1* ','C8  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,7,'C6  ','N3  ','N1  ','C4  ',2,'C8  ','C1* ',
     1  'O6  ','N7  ','T','S:W','W:S',1)

c       QUO-QUO (WST) HB betn. N1G-N3G & O6G-C1*G
c       call findpair(i,7,'C4  ','N1  ','N3  ','C6  ',7,'C7  ','O6  ',
c    1  'C1* ','C8  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,7,'C6  ','N3  ','N1  ','C4  ',7,'C8  ','C1* ',
     1  'O6  ','C7  ','T','S:W','W:S',1)

c       GUA-N6G (WST) HB betn. N1G-N3G & O6G-C1*G
c       call findpair(i,2,'C4  ','N1  ','N3  ','C6  ',5,'N7  ','O6  ',
c    1  'C1* ','C8  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,5,'C6  ','N3  ','N1  ','C4  ',2,'C8  ','C1* ',
     1  'O6  ','N7  ','T','S:W','W:S',1)

c       QUO-N6G (WST) HB betn. N1G-N3G & O6G-C1*G
c       call findpair(i,7,'C4  ','N1  ','N3  ','C6  ',5,'C7  ','O6  ',
c    1  'C1* ','C8  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,5,'C6  ','N3  ','N1  ','C4  ',7,'C8  ','C1* ',
     1  'O6  ','C7  ','T','S:W','W:S',1)

c       GUA-GUA (HSC) HB betn. N7G-C1*G & O6G-N2G 
c       call findpair(i,2,'N1  ','O6  ','N2  ','N1  ',2,'N9  ','N7  ',
c    1  'C1* ','C8  ','C','H:S','S:H',1)    Sukanya
c       call findpair(i,2,'N1  ','N2  ','O6  ','N1  ',2,'C8  ','C1* ',
c    1  'N7  ','N9  ','C','S:H','H:S',1) removed by Sukanya 
c        1 structure in 1VQO, buckle= -70

c       GUA-QUO (HSC) HB betn. N7G-C1*G & O6G-N2G 
c       call findpair(i,2,'N1  ','O6  ','N2  ','N1  ',7,'N9  ','N7  ',
c    1  'C1* ','C8  ','C','H:S','S:H',1)    Sukanya
        call findpair(i,7,'N1  ','N2  ','O6  ','N1  ',2,'C8  ','C1* ',
     1  'N7  ','N9  ','C','S:H','H:S',1)

c       GUA-N6G (HSC) HB betn. N7G-C1*G & O6G-N2G 
c       call findpair(i,2,'N1  ','O6  ','N2  ','N1  ',5,'N9  ','N7  ',
c    1  'C1* ','C8  ','C','H:S','S:H',1)    Sukanya
        call findpair(i,5,'N1  ','N2  ','O6  ','N1  ',2,'C8  ','C1* ',
     1  'N7  ','N9  ','C','S:H','H:S',1)

c       GUA-GUA (HST) HB betn. N7G-N2G & O6G-C1*G
c       call findpair(i,2,'C4  ','N7  ','N2  ','N1  ',2,'N3  ','O6  ',
c    1  'C1* ','C8  ','T','H:S','S:H',1)    Sukanya
        call findpair(i,2,'N1  ','N2  ','N7  ','C4  ',2,'C8  ','C1* ',
     1  'O6  ','N3  ','T','S:H','H:S',1)

c       GUA-GUA (SSC) HB betn. N2G-N3G & N3G-C1*G
c       call findpair(i,2,'N1  ','N2  ','N3  ','C5  ',2,'C6  ','N3  ',
c    1  'C1* ','C8  ','C','S:S','S:S',1)    Sukanya
        call findpair(i,2,'C5  ','N3  ','N2  ','N1  ',2,'C8  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c*******************************GUA-CYT********************************

c       GUA-CYT (WSC) HB betn. N1G-O2C & N2G-C1*C
        call findpair(i,3,'N3  ','O2  ','N1  ','C5  ',2,'C6  ','C1* ',
     1  'N2  ','N3  ','C','S:W','W:S',1)

c       QUO-CYT (WSC) HB betn. N1G-O2C & N2G-C1*C
        call findpair(i,3,'N3  ','O2  ','N1  ','C5  ',7,'C6  ','C1* ',
     1  'N2  ','N3  ','C','S:W','W:S',1)

c       N6G-CYT (WSC) HB betn. N1G-O2C & N2G-C1*C
        call findpair(i,3,'N3  ','O2  ','N1  ','C5  ',5,'C6  ','C1* ',
     1  'N2  ','N3  ','C','S:W','W:S',1)

c       GUA-CYT (WST) 1st HB betn. N2G-O2C & N1G-C1*C
        call findpair(i,3,'N3  ','O2  ','N2  ','N3  ',2,'C6  ','C1* ',
     1  'N1  ','C5  ','T','S:W','W:S',1)

c       QUO-CYT (WST) 1st HB betn. N2G-O2C & N1G-C1*C
        call findpair(i,3,'N3  ','O2  ','N2  ','N3  ',7,'C6  ','C1* ',
     1  'N1  ','C5  ','T','S:W','W:S',1)

c       N6G-CYT (WST) HB betn. N2G-O2C & N1G-C1*C
        call findpair(i,3,'N3  ','O2  ','N2  ','N3  ',5,'C6  ','C1* ',
     1  'N1  ','C5  ','T','S:W','W:S',1)

c       GUA-CYT (WST) 2nd HB betn. N1G-O2C & O6G-C1*C
c       call findpair(i,2,'C4  ','N1  ','O2  ','C5  ',3,'N7  ','O6  ',
c    1  'C1* ','C6  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,3,'C5  ','O2  ','N1  ','C4  ',2,'C6  ','C1* ',
     1  'O6  ','N7  ','T','S:W','W:S',1)

c       QUO-CYT (WST) 2nd HB betn. N1G-O2C & O6G-C1*C
c       call findpair(i,7,'C4  ','N1  ','O2  ','C5  ',3,'C7  ','O6  ',
c    1  'C1* ','C6  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,3,'C5  ','O2  ','N1  ','C4  ',7,'C6  ','C1* ',
     1  'O6  ','C7  ','T','S:W','W:S',1)

c       GUA-CYT (SWC) HB betn. N3G-N4C & C1*G-N3C
        call findpair(i,2,'N7  ','N3  ','N4  ','C5  ',3,'C8  ','C1* ',
     1  'N3  ','C6  ','C','S:W','W:S',1)

c       QUO-CYT (SWC) HB betn. N3G-N4C & C1*G-N3C
        call findpair(i,7,'C7  ','N3  ','N4  ','C5  ',3,'C8  ','C1* ',
     1  'N3  ','C6  ','C','S:W','W:S',1)

c       GUA-CYT (SSC) HB betn. N2G-O2C & N3G-C1*C
        call findpair(i,3,'N3  ','O2  ','N2  ','N1  ',2,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c       QUO-CYT (SSC) HB betn. N2G-O2C & N3G-C1*C
        call findpair(i,3,'N3  ','O2  ','N2  ','N1  ',7,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c       N6G-CYT (SSC) HB betn. N2G-O2C & N3G-C1*C
        call findpair(i,3,'N3  ','O2  ','N2  ','N1  ',5,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c       GUA-CYT (SST) HB betn. N2G-C1*C & C1*G-O2C
        call findpair(i,3,'C5  ','O2  ','C1* ','C8  ',2,'C6  ','C1* ',
     1  'N2  ','N1  ','T','S:S','S:S',1)

c*******************************GUA-URA********************************

c       GUA-URA (WSC) HB betn. N1G-O2U & N2G-C1*U
        call findpair(i,4,'N3  ','O2  ','N1  ','C5  ',2,'C6  ','C1* ',
     1  'N2  ','N3  ','C','S:W','W:S',1)

c       QUO-URA (WSC) HB betn. N1G-O2U & N2G-C1*U
        call findpair(i,4,'N3  ','O2  ','N1  ','C5  ',7,'C6  ','C1* ',
     1  'N2  ','N3  ','C','S:W','W:S',1)

c       N6G-URA (WSC) HB betn. N1G-O2U & N2G-C1*U
        call findpair(i,4,'N3  ','O2  ','N1  ','C5  ',5,'C6  ','C1* ',
     1  'N2  ','N3  ','C','S:W','W:S',1)

c       GUA-URA (WST) HB betn. N2G-O2U & N1G-C1*U
        call findpair(i,4,'N3  ','O2  ','N2  ','N3  ',2,'C6  ','C1* ',
     1  'N1  ','C5  ','T','S:W','W:S',1)

c       QUO-URA (WST) HB betn. N2G-O2U & N1G-C1*U
        call findpair(i,4,'N3  ','O2  ','N2  ','N3  ',7,'C6  ','C1* ',
     1  'N1  ','C5  ','T','S:W','W:S',1)

c       N6G-URA (WST) HB betn. N2G-O2U & N1G-C1*U
        call findpair(i,4,'N3  ','O2  ','N2  ','N3  ',5,'C6  ','C1* ',
     1  'N1  ','C5  ','T','S:W','W:S',1)

c       GUA-PSU (SHC) HB betn. N2G-O2U & C1*G-N1U
        call findpair(i,2,'N1  ','N2  ','O2  ','N3  ',6,'C8  ','C1* ',
     1  'N1  ','C5  ','C','S:H','H:S',1)

c       QUO-PSU (SHC) HB betn. N2G-O2U & C1*G-N1U
        call findpair(i,7,'N1  ','N2  ','O2  ','N3  ',6,'C8  ','C1* ',
     1  'N1  ','C5  ','C','S:H','H:S',1)

c       N6G-PSU (SHC) HB betn. N2G-O2U & C1*G-N1U
        call findpair(i,5,'N1  ','N2  ','O2  ','N3  ',6,'C8  ','C1* ',
     1  'N1  ','C5  ','C','S:H','H:S',1)

c       GUA-URA (SSC) HB betn. N2G-O2U & N3G-C1*U
        call findpair(i,4,'N3  ','O2  ','N2  ','N1  ',2,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c       QUO-URA (SSC) HB betn. N2G-O2U & N3G-C1*U
        call findpair(i,4,'N3  ','O2  ','N2  ','N1  ',7,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c       N6G-URA (SSC) HB betn. N2G-O2U & N3G-C1*U
        call findpair(i,4,'N3  ','O2  ','N2  ','N1  ',5,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c       GUA-URA (SST) HB betn. N2G-C1*U & C1*G-O2U (to be ignored,
c due to two H-bonds involving sugar.   bhatta March 2020
c        call findpair(i,4,'C5  ','O2  ','C1* ','C8  ',2,'C6  ','C1* ',
c     1  'N2  ','N1  ','T','S:S','S:S',1)

c       GUA-PSU (SSC) HB betn. N2G-O4U & N3G-C1*U
c       call findpair(i,2,'N1  ','N2  ','O4  ','C2  ',6,'C6  ','N3  ',
c    1  'C1* ','C6  ','C','S:S','S:S',1)    Sukanya
        call findpair(i,6,'C2  ','O4  ','N2  ','N1  ',2,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c       QUO-PSU (SSC) HB betn. N2G-O4U & N3G-C1*U
c       call findpair(i,7,'N1  ','N2  ','O4  ','C2  ',6,'C6  ','N3  ',
c    1  'C1* ','C6  ','C','S:S','S:S',1)    Sukanya
        call findpair(i,6,'C2  ','O4  ','N2  ','N1  ',7,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c       N6G-PSU (SSC) HB betn. N2G-O4U & N3G-C1*U
c       call findpair(i,5,'N1  ','N2  ','O4  ','C2  ',6,'C6  ','N3  ',
c    1  'C1* ','C6  ','C','S:S','S:S',1)    Sukanya
        call findpair(i,6,'C2  ','O4  ','N2  ','N1  ',5,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:S','S:S',1)

c*******************************CYT-CYT********************************

c       CYT-CYT (SWC) HB betn. O2C-N4C & C1*C-N3C
        call findpair(i,3,'N3  ','O2  ','N4  ','C5  ',3,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:W','W:S',1)

c       CYT-CYT (S+T) HB betn. O2C-N3C & C1*C-N4C
        call findpair(i,3,'N3  ','O2  ','N3  ','C6  ',3,'C6  ','C1* ',
     1  'N4  ','C5  ','T','S:+','+:S',1)

c       CYT-CYT (SST) HB betn. O2C-C1*C & C1*C-O2C
        call findpair(i,3,'C5  ','O2  ','C1* ','C6  ',3,'C6  ','C1* ',
     1  'O2  ','C5  ','T','S:S','S:S',1)

c*******************************CYT-URA********************************

c       CYT-URA (WSC) HB betn. N4C-O2U & N3C-C1*U
        call findpair(i,4,'N3  ','O2  ','N4  ','C5  ',3,'C6  ','C1* ',
     1  'N3  ','C6  ','C','S:W','W:S',1)

c       CYT-PSU (WSC) HB betn. N4C-O4U & O2C-C1*U
c       call findpair(i,3,'C5  ','N4  ','O4  ','N1  ',6,'N1  ','O2  ',
c    1  'C1* ','C6  ','C','W:S','S:W',1)    Sukanya
        call findpair(i,6,'N1  ','O4  ','N4  ','C5  ',3,'C6  ','C1* ',
     1  'O2  ','N1  ','C','S:W','W:S',1)

c       CYT-URA (+ST) HB betn. N3C-O2U & N4C-C1*U
        call findpair(i,4,'N3  ','O2  ','N3  ','C6  ',3,'C6  ','C1* ',
     1  'N4  ','C5  ','T','S:+','+:S',1)

c       CYT-URA (SWC) HB betn. O2C-N3U & C1*C-O2U
        call findpair(i,3,'C5  ','O2  ','N3  ','C6  ',4,'C6  ','C1* ',
     1  'O2  ','C1* ','C','S:W','W:S',1)

c       CYT-URA (SWT) HB betn. O2C-N3U & C1*C-O4U
        call findpair(i,3,'C5  ','O2  ','N3  ','C6  ',4,'C6  ','C1* ',
     1  'O4  ','C5  ','T','S:W','W:S',1)

c       CYT-PSU (SHT) HB betn. O2C-N1U & C1*C-O2U
        call findpair(i,3,'C5  ','O2  ','N1  ','C4  ',6,'C5  ','N4  ',
     1  'O2  ','N3  ','T','S:H','H:S',1)

c       CYT-URA (SST) HB betn. O2C-C1*U & C1*C-O2U
        call findpair(i,3,'C5  ','O2  ','C1* ','C6  ',4,'C6  ','C1* ',
     1  'O2  ','C5  ','T','S:S','S:S',1)

c*******************************URA-URA********************************

c       URA-URA (WSC) HB betn. N3U-O2U & O2U-C1*U
c       call findpair(i,4,'C6  ','N3  ','O2  ','C5  ',4,'N1  ','O2  ',
c    1  'C1* ','C6  ','C','W:S','S:W',1)    Sukanya
        call findpair(i,4,'C5  ','O2  ','N3  ','C6  ',4,'C6  ','C1* ',
     1  'O2  ','N1  ','C','S:W','W:S',1)

c       URA-PSU (WSC) HB betn. N3U-O4U & O2U-C1*U
c       call findpair(i,4,'C6  ','N3  ','O4  ','C2  ',6,'N1  ','O2  ',
c    1  'C1* ','C6  ','C','W:S','S:W',1)    Sukanya
        call findpair(i,6,'C2  ','O4  ','N3  ','C6  ',4,'C6  ','C1* ',
     1  'O2  ','N1  ','C','S:W','W:S',1)

c       URA-URA (WST) HB betn. N3U-O2U & O4U-C1*U
c       call findpair(i,4,'C6  ','N3  ','O2  ','C5  ',4,'C5  ','O4  ',
c    1  'C1* ','C6  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,4,'C5  ','O2  ','N3  ','C6  ',4,'C6  ','C1* ',
     1  'O4  ','C5  ','T','S:W','W:S',1)

c       URA-PSU (WST) HB betn. N3U-O4U & O4U-C1*U
c       call findpair(i,4,'C6  ','N3  ','O4  ','C2  ',6,'C5  ','O4  ',
c    1  'C1* ','C6  ','T','W:S','S:W',1)    Sukanya
        call findpair(i,6,'C2  ','O4  ','N3  ','C6  ',4,'C6  ','C1* ',
     1  'O4  ','C5  ','T','S:W','W:S',1)

c       URA-PSU (SHT) HB betn. O2U-N1U & C1*U-O2U
        call findpair(i,4,'C5  ','O2  ','N1  ','C4  ',6,'C6  ','C1* ',
     1  'O2  ','N3  ','T','S:H','H:S',1)

c       URA-URA (SST) HB betn. O2U-C1*U & C1*U-O2U
        call findpair(i,4,'C5  ','O2  ','C1* ','C6  ',4,'C6  ','C1* ',
     1  'O2  ','C5  ','T','S:S','S:S',1)

c**********************************************************************

	endif

c        write(*,*)'H-bond found',i
        enddo

c       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
c       ----------------------------------------------------------###
c       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         do i = 1,nres
          do jj=1,10
            tpod(jj)=' '
          end do
          do j = 1,base(1,i)
	   localinfo(j)=base(j+1,i)
           energylc(j)=energy(j,i)
	   bpinfo(j)=feature(j,i)
	   tpinfo(j)=type(j,i)
           anglinfo(j)=yangle(j,i)
	   bptype(j)=bpinfo(j)
c 
c Converting different types of protonated basepair edges into their equivalent
c ones
c
	   if(bptype(j)(1:1).eq.'+') bptype(j)(1:1)='W'
	   if(bptype(j)(1:1).eq.'w') bptype(j)(1:1)='W'
	   if(bptype(j)(3:3).eq.'+') bptype(j)(3:3)='W'
	   if(bptype(j)(3:3).eq.'w') bptype(j)(3:3)='W'
	   if(bptype(j)(3:3).eq.'s') bptype(j)(3:3)='S'
	   if(bptype(j)(3:3).eq.'z') bptype(j)(3:3)='S'
	   if(bptype(j)(3:3).eq.'h') bptype(j)(3:3)='H'
	   if(bptype(j)(1:1).eq.'h') bptype(j)(1:1)='H'
	   if(bptype(j)(1:1).eq.'s') bptype(j)(1:1)='S'
	   if(bptype(j)(1:1).eq.'z') bptype(j)(1:1)='S'
c********************************************************************
c   TO KEEP ONLY WATSON CRICK
c********************************************************************
           if(bptype(j).ne.'W:W')then
               tpod(j)='O'
           elseif(bptype(j).eq.'W:W')then
             if(bpinfo(j).eq.'W:W')then
                tpod(j)='W'
             else
                tpod(j)='N'
             endif
           endif
	  enddo
          nj=0
          nc=0
          do j=1,base(1,i)
            if(base(1,i).le.1)then
              exit
            else
            if(tpod(j).eq.'N')then
              do jl=1,10
               if(tpod(jl).eq.'W') nc=nc+1 
              end do 
              if(nc.gt.0)then 
                nj=nj+1
                if(tpod(j+1).ne.' ')then
                nn=j
c
c   Keeping only the main H-bonding edge if both canonical form and protonated
c   forms of a base are found to be Hbonding to another base with its same
c   edge.
c
                do jk=nn,10
                  bpinfo(jk)=bpinfo(jk+1)
                  tpinfo(jk)=tpinfo(jk+1)
                  energylc(jk)=energylc(jk+1)
                  localinfo(jk)=localinfo(jk+1)
                  bptype(jk)=bptype(jk+1)
                  anglinfo(jk)=anglinfo(jk+1)
                end do
                elseif(tpod(j+1).eq.' ')then
                  bpinfo(j)='   '
                  tpinfo(j)=' '
                  energylc(j)=0.0            
                  localinfo(j)=0               
                  bptype(j)='   '        
                  anglinfo(j)=0.0             
                endif
               endif
              endif
            endif
          end do 
          base(1,i)=base(1,i)-nj             
c          if(base(1,i).eq.0)then
c	  elseif(base(1,i).ge.1) then
	  if(base(1,i).ge.1) then
	      nopair=base(1,i)
              do j = 1,nopair
	        L=J+1
                do while (L.le.base(1,i))
c
c   Sorting to bring the energetically best base pair as the principal one by
c   pushing the other ones as secondary type (BF or TP)
c
	          if(energylc(j).gt.energylc(L))then
                     bsinfo=localinfo(L) 
                     localinfo(L)=localinfo(j)
                     localinfo(j)=bsinfo
                     bs=bpinfo(L)
                     bpinfo(L)=bpinfo(j)
                     bpinfo(j)=bs
                     bs=bptype(L)
                     bptype(L)=bptype(j)
                     bptype(j)=bs
                     tp=tpinfo(L)
                     tpinfo(L)=tpinfo(j)
                     tpinfo(j)=tp 
                     varenergy=energylc(L)
                     energylc(L)=energylc(j)
                     energylc(j)=varenergy
                     varangl=anglinfo(L)
                     anglinfo(L)=anglinfo(j)
                     anglinfo(j)=varangl
                  endif 
	          L=L+1
	        enddo
	      enddo
c
c   #######  module for bifurcated hydrogen bonded basepairs ######
c
            do j = 1,nopair
	      L=J+1
              do while (L.le.base(1,i))
               if(localinfo(j).eq.localinfo(l).or.  
     1            bptype(j)(1:1).eq.bptype(l)(1:1)) then
	          do k = l,base(1,i)-1
	            localinfo(k)=localinfo(k+1)
	            bpinfo(k)=bpinfo(k+1)
	            bptype(k)=bptype(k+1)
	            tpinfo(k)=tpinfo(k+1)
                    energylc(k)=energylc(k+1)
                    anglinfo(k)=anglinfo(k+1)
                  enddo
	          base(1,i)=base(1,i)-1
	       else
	          L=L+1
	       endif
              enddo
            enddo
          endif
	  do j=1,base(1,i)
            l=j+1
c           if(l.le.base(1,i))then
              base(l,i)=localinfo(j)
c           endif
            energy(j,i)=energylc(j)
            feature(j,i)=bpinfo(j)
            type(j,i)=tpinfo(j)
            yangle(j,i)=anglinfo(j)
          end do
         enddo
c----------------------------------------------------------------
c        TO sort into BF and TRP
c---------------------------------------------------------------
	do i=1,nres
         do j=1,base(1,i)
            l=j+1
            do k=1,nres
              if(k.ne.i)then
                ipair=base(l,i)
                if(ipair.eq.k)then
                  do m=1,base(1,k)
                     n=m+1
                     kpair=base(n,k)
                     if(kpair.eq.i)then
                        ncount=1
                        exit
                     else
                       ncount=0
                     endif
                  end do
                  if(ncount.eq.0)then
                    num=base(1,k)
                    base(num+2,k)=i
                    base(1,k)=num+1
                    energy(num+1,k)=energy(j,i)
                    feature(num+1,k)(1:1)=feature(j,i)(3:3)
                    feature(num+1,k)(3:3)=feature(j,i)(1:1)
                    feature(num+1,k)(2:2)=feature(j,i)(2:2)
                    type(num+1,k)=type(j,i)
                    yangle(num+1,k)=yangle(j,i)
                  endif
                endif
              endif
            end do
         end do          
        end do
c        do i=1,nres
c          if(base(1,i).gt.1)then
c            do j=1,base(1,i)
c              L=J+1
c              if(L.le.base(1,i))then
c                if(energy(j,i).gt.energy(l,i))then
c                   bsinfo=base(l+1,i) 
c                   base(l+1,i)=base(l,i)   
c                   base(l,i)=bsinfo
c                   bs=feature(l,i)
c                   feature(l,i)=feature(j,i)
c                   feature(j,i)=bs
c                   tp=type(l,i)
c                   type(l,i)=type(j,i)
c                   type(j,i)=tp 
c                   varenergy=energy(l,i)
c                   energy(l,i)=energy(j,i)
c                   energy(j,i)=varenergy
c                   varangl=yangle(L,i)
c                   yangle(L,i)=yangle(j,i)
c                   yangle(j,i)=varangl
c                endif         
c              endif
c            end do
c          endif
        do i=1,nres
         if(i.gt.1) then
          if(base(1,i).gt.1) then
             do j=2,base(1,i)
               L=j+1
c
c   Instead of original energy based sorting to detect the Primary Paired
c   base (BP) now it would try to sort depending on the previous basepair''s 
c   residue number.
c
                if(((base(2,i-1).eq.(base(l,i)+1)).or.(base(2,i-1).eq.
     1              (base(l,i)-1))).and.base(2,i-1).ne.0) then
c       write(*,*) 'IF CONDITION: I, J',i,j,base(1,i)
                   bsinfo=base(l,i)
                   base(l,i)=base(j,i)
                   base(j,i)=bsinfo
                   bs=feature(l-1,i)
                   feature(l-1,i)=feature(j-1,i)
                   feature(j-1,i)=bs
                   tp=type(l-1,i)
                   type(l-1,i)=type(j-1,i)
                   type(j-1,i)=tp
                   varenergy=energy(l-1,i)
                   energy(l-1,i)=energy(j-1,i)
                   energy(j-1,i)=varenergy
                   varangl=yangle(l-1,i)
                   yangle(l-1,i)=yangle(j-1,i)
                   yangle(j-1,i)=varangl
                else if(((base(2,i+1).eq.(base(l,i)+1)).or.(base(2,i+1)
     1     .eq.(base(l,i)-1))).and.base(2,i+1).ne.0.and.i.lt.nres.and.
     2     base(1,i-1).eq.1) then
c       write(*,*) 'ELSEIF CONDITION: I, J',i,j,base(1,i)
c
c   After sorting Triplet steps based on the previous basepair''s residue number
c   it would now sort according the the next residue no.
c
                   bsinfo=base(l,i)
                   base(l,i)=base(j,i)
                   base(j,i)=bsinfo
                   bs=feature(l-1,i)
                   feature(l-1,i)=feature(j-1,i)
                   feature(j-1,i)=bs
                   tp=type(l-1,i)
                   type(l-1,i)=type(j-1,i)
                   type(j-1,i)=tp
                   varenergy=energy(l-1,i)
                   energy(l-1,i)=energy(j-1,i)
                   energy(j-1,i)=varenergy
                   varangl=yangle(l-1,i)
                   yangle(l-1,i)=yangle(j-1,i)
                   yangle(j-1,i)=varangl
                 endif
             enddo
           endif
         endif

          do j=1,10       
             bptp(j)="  "
          end do 
          do j=1,base(1,i)
	     bptype(j)=feature(j,i)
	     if(bptype(j)(1:1).eq.'+') bptype(j)(1:1)='W'
	     if(bptype(j)(1:1).eq.'w') bptype(j)(1:1)='W'
	     if(bptype(j)(3:3).eq.'+') bptype(j)(3:3)='W'
	     if(bptype(j)(3:3).eq.'w') bptype(j)(3:3)='W'
	     if(bptype(j)(3:3).eq.'s') bptype(j)(3:3)='S'
	     if(bptype(j)(3:3).eq.'z') bptype(j)(3:3)='S'
	     if(bptype(j)(3:3).eq.'h') bptype(j)(3:3)='H'
	     if(bptype(j)(1:1).eq.'h') bptype(j)(1:1)='H'
	     if(bptype(j)(1:1).eq.'s') bptype(j)(1:1)='S'
	     if(bptype(j)(1:1).eq.'z') bptype(j)(1:1)='S'
          end do
          if(base(1,i).ge.1)then
            do j=1,base(1,i)
              bptp(1)="BP"
              L=J+1
              if(bptype(1)(1:1).eq.bptype(l)(1:1))then
                bptp(l)="BF"
              else
                bptp(l)="TP"
              endif
            end do
          endif
          if(base(1,i).eq.0)then
             nchk=0
          else
             do l=1,base(1,i)
               if(bptp(l).eq."BF")then
                 nchk=1
                 newj=l
                 exit
               else
                 nchk=0
               endif
             end do
          endif
c           write(*,*) 'Final',i,base(2,i),feature(1,i),type(1,i),
c     1 energy(1,i),yangle(1,i),' BASE(1,i)',base(1,i)

c	if(oligo.eq.1.and.(base(2,i).eq.(i-1).and.base(2,i-1).eq.i)) then
c	   ierror=1
c	   return
c	endif  		! Commented out module, bhatta March 17, 2021

!         write(*,331) i,prd(i)
331    format(I5,I5)
C        write(6,*) i,base(1,i),base(2,i),base(3,i)
          do kpr=1,40
            prnvar(kpr)=''
          enddo
c        write(6,*) i,prd(i)
          write(prnvar(1),*) i
          write(prnvar(2),*) prd(i)
          write(prnvar(3),*) resd(ists(i))
          write(prnvar(4),*) allins(ists(i))
          write(prnvar(5),*) pcd(i)
          prntins=allins(ists(j))
c        write(6,*) i,'PRNTINS ',prntins,' ',prnvar(5)
c        write(6,*) 'base(1,i)=',base(1,i)
          do j=1,base(1,i)
            write(prnvar((j-1)*9+6),*) base(j+1,i)
            write(prnvar((j-1)*9+7),*) prd(base(j+1,i))
            write(prnvar((j-1)*9+8),*) resd(ists(base(j+1,i)))
            write(prnvar((j-1)*9+9),*) allins(ists(base(j+1,i)))
            write(prnvar((j-1)*9+10),*) pcd(base(j+1,i))
            write(prnvar((j-1)*9+11),*) feature(j,i)
            write(prnvar((j-1)*9+12),*) type(j,i)
            write(prnvar((j-1)*9+13),*) bptp(j)
            write(prnvar((j-1)*9+14),'(f10.2)') energy(j,i)
          enddo
          do kpr=1,40
            prnvar(kpr)=adjustl(prnvar(kpr))
            nnfpr(kpr)=index(prnvar(kpr),' ')
          enddo
c        write(53,*)nmpass(1:4),' ',(prnvar(kpr)(1:nnfpr(kpr)),kpr=1,36),
c     1   prntins

	if(nocsv.eq.0) then
        write(53,*)(adjustl(nmpass(1:4))),' ',(prnvar(kpr)(1:nnfpr(kpr))
     1   ,kpr=1,3),(prnvar(kpr)(1:nnfpr(kpr)), kpr=4,36)
        nnf=index(line,'  ')
c        write(53,5) line(2:nnf-1)
	endif

c        write(6,*) i,prd(i),resd(ists(i)),allins(ists(i)),pcd(i)
	if(nevalue.eq.0) then
 	   write(52,116) i,prd(i),resd(ists(i)),allins(ists(i)),pcd(i),
     1 (base(j+1,i),prd(base(j+1,i)),
     2  resd(ists(base(j+1,i))),allins(ists(base(j+1,i))),
     2 pcd(base(j+1,i)),feature(j,i),type(j,i),bptp(j),energy(j,i),
     3 j=1,base(1,i))

	else
	  do j=1,base(1,i)
	    k=base(j+1,i)
	    if(base(j+1,i).gt.0) val2pr(j)=calcc1dis(i,k)
	  enddo
	  
 	   write(52,119) i,prd(i),resd(ists(i)),allins(ists(i)),pcd(i),
     1 (base(j+1,i),prd(base(j+1,i)),
     2  resd(ists(base(j+1,i))),allins(ists(base(j+1,i))),
     2 pcd(base(j+1,i)),feature(j,i),type(j,i),bptp(j),
     3 val2pr(j), j=1,base(1,i))

	endif
        line=''

            if(feature(1,i).eq.'W:W') strseq(i)='W'
            if(feature(1,i).ne.'W:W') strseq(i)='N'
            if(base(1,i).eq.0) strseq(i)='C'
            if(base(1,i).gt.1) strseq(i)='T'
            if(nchk.eq.0)then
c 
c  This statement stops writing mirror image section of nup file appearing
c  due to antiparallel double helical nature of regular oligonucleotides.
c
	    if(nonup.eq.0) then
              write(62,165) i,(base(j+1,i),feature(j,i),type(j,i),
     1        bptp(j),j=1,base(1,i))
	    endif
            elseif(nchk.eq.1)then
              do k=1,4
                locp(k)=0        
                locf(k)="   "
                loct(k)=" "
                locb(k)="  "
              end do
              do k=1,base(1,i)
                if(k.lt.newj)then
                  locp(k)=base(k+1,i)
                  locf(k)=feature(k,i)
                  loct(k)=type(k,i)
                  locb(k)=bptp(k)
                elseif(k.eq.newj)then
                   locp(4)=base(k+1,i)
                   locf(4)=feature(k,i)
                   loct(4)=type(k,i)
                   locb(4)=bptp(k)
                elseif(k.gt.newj)then
                   locp(k-1)=base(k+1,i)
                   locf(k-1)=feature(k,i)
                   loct(k-1)=type(k,i)
                   locb(k-1)=bptp(k)
                endif
              end do
	      if(nonup.eq.0) then
              write(62,165) i,(locp(k),locf(k),loct(k),locb(k),k=1,4)
	      endif
           endif  
          end do 
        close(unit=4)
        close(unit=9)
        close(unit=52)
        close(unit=53)
        close(unit=62)
        close(unit=77)
!        close(unit=79)  ! done by parthajit roy. I don't think bhatta
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TO FIND OUT THE LOOP & BULDGE REGIONS
! Debasish Mukherjee 8th April 2014
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!        icnt=0
!        j=0
!        ik=0
!        open(unit=10,file=nmpass)
!        do i=1,10000
!        read(10,5,end=100) aline
!         if(aline(1:1) .eq. ' ')then
!         read(aline,101) iaa(i),iba(i),aba(i),aca(i),ica(i),ida(i),
!     1 baa(i),bca(i)
!           if(ica(i) .eq. 0) then
!           icnt=icnt+1
!           else
!           iloop(i)=icnt
!           icnt=0
!           endif
!         else
!         ik=ik+1
!         endif
!        enddo
100     continue
101   format(I5,I5,1x,A3,1x,a1,I5,I5,1x,A3,1x,a1)
!        do j=1,i
!          if(iloop(j) .gt. 0) then
!            if(((iba(j)-iloop(j)-1) .eq. ida(j)) .and.(aca(j) .eq. 
!     1 bca(j)))then
!            do k=1,iloop(j)
!            l=j-iloop(j)-ik-1+k
!            strseq(l)='C'
!            enddo
!            else
!            do k=1,iloop(j)
!            l=j-iloop(j)-ik-1+k
!	    strseq(l)='B'
!	    enddo
!	    endif
!	  endif
!	enddo
	!======================================================================
c	  write(22,168) (strseq(i),i=1,nres)
c        write(*,*) 'I and NRES',i,nres
c4       format(13X,A4,A3,1X,A1,I4,a1,3X,3F8.3,2x,f4.2,f6.2)
c	if(i.eq.nres+1) then
        call rnahelix(nmpass,filenm,nnf,modelno,nocsv,nodat,nodbn,nohlx)
c        write(6,*) 'Calculation of a MODEL over'
!        enddo                                   ! end of kmodel loop
c        endif
	return
211	continue
	write(*,*) 'Error opening input file:',filenm
	ierror=1
        return
256     continue
        write(6,*)'Error opening database files, AdeVariants.name, etc.'
        write(6,*) 'Please make sure to keep these files in a directory 
     1and pass that absolute path through NUCLEIC_ACID_DIR environment
     2variable using Unix Shell command.'
	ierror=1
	return
5       format(A80)
255       format('#',A80)
4       format(12X,A4,1x,A3,1X,A1,I4,a1,3X,3F8.3,2x,f4.2,f6.2)
6     format('#HEADER   Base Pairing/Triplet Scheme for Structure',
     1 ' File',a20)
13      format(I5,i5,1x,A1,I5)
16      format(I5,I5,1x,A1,4(I5,I5,1x,A1,2x,A3,1X,A1,f8.5))
77    format('#====================================================='
     1  ,'==========================')
60    format('#Column starting from > mark represent the following:')
61    format('#>Serial number of nucleotides')
62    format('#',6('-'),'>Residue number as in PDB')
63    format('#',16('-'),'>Residue name')
160   format('#====='/'# Generated by BPFIND release June 2020, version
     1 2.1.17'/'#=====')
163   format('#',18('-'),'>PDB_ins_code. Useful in tRNA with long',
     1' variable loop')
C This is useful for ',
C     1'unique identification of additional residues, 47 in tRNA variable
C     2 loop, for example')
64    format('#',20('-'),'>Chain ID')
65    format('#',22('-'),'>Serial number of the paired base')
66    format('#',30('-'),'>Residue number of the paired base as in PDB')
67    format('#',39('-'),'>Residue name of the paired base')
164   format('#',41('-'),'>PDB_ins_code of the paired residue')
68    format('#',43('-'),'>Chain ID of the paired base')
69    format('#',47('-'),'>Base pair type')
70    format('#',53('-'),'>Base pair indicator')
71    format('#',55('-'),'>E-value indicating base pair deformation')
72    format('#',62('-'),'>Equivalent information for other pairs, such
     1as in base triples') 
73    format('#HEADER   Program BPFIND run with the following options:')
74    format('#HEADER   Distance cutoff: ',f3.1,'   Angle cutoff: ',
     1         f5.1,'   E-value cutoff: ',f3.1)
75    format('#HEADER   Results for NMR model no.',i6)
275     format('#HEADER   Consideration of CH...O H-Bond: ',a3)    
276     format('#HEADER   Consideration of Sugar mediated H-bondB: ',a3)
277     format('#HEADER   Choise of HETATM entries: ',a3)   ! ans(hetatm+1)
278     format('#HEADER   Option as Oligonucleotide: ',a3)  ! ans(oligo+1)
375     format('#HEADER   HT option was selected to consider non-natural
     1 nucleotides also')
376     format('#HEADER   CH option was selected to avoid identification 
     1 of base pairs involving C-H...O interaction')
377     format('#HEADER   SG option was selected to avoid identification
     1 of base pairs involving 2''OH group of sugar')
79      format(8x,i8)
116   format(1x,I5,1x,I7,1x,A3,1x,a1,1x,A4,4(1x,I5,1x,I5,1x,
     1   a3,1x,a1,1x,A4,1x,A3,A1,1x,a2,f5.2))

119   format(1x,I5,1x,I7,1x,A3,1x,a1,1x,A4,4(1x,I5,1x,I5,1x,
     1   a3,1x,a1,1x,A4,1x,A3,A1,1x,a2,f6.2))

165     format(I5,1x,4(I5,1x,A3,1X,A1,1x,a2))
167     format(a5) 
168   format(70a1)

210     write(6,9999)
	call showoptn(ierror)
c	stop 1
        end

C      ----------------------------------------------------------######

C      ----------------------------------------------------------######
       subroutine findhb3(ires,ib1,ib2,atp11,at11,at21,atp21,atp12,
     1   at12,at22,atp22,atp13,at13,at23,atp23,config,
     2   hbftr,hbftrbck)
       EXTERNAL VERIFY
       integer ib1,ib2,ATGC,KS,NN1,NN2,ires,atmN,base,NN3,NN4,KS11,KS12
     1          ,KS21,KS22,KS31,KS32,ncback
       real   distcut,anglcut,energy,calenrg12,anglconv,yangle,
     1 calenrg13,calenrg,rangct,occ, invdegree
       character*4 at11,atp11,at21,at12,atp12,at22,named,atp21,
     1              atp22,at13,atp13,at23,atp23
       character*3 hbftr,feature,hbftrbck,resid,resd
       character*1 config,type,allins
       character*4 chaind 
       COMMON /HBQUA/BASE(21,20000),energy(20,20000),
     1                    yangle(20,20000)
       COMMON /HBQTYP/TYPE(20,20000),FEATURE(20,20000)
       COMMON /NAMES/NAMED(900000),RESD(900000)
       COMMON /CHAINS/CHAIND(900000)
       COMMON /ALPHA/XD(900000),YD(900000),ZD(900000)
       COMMON /BETA/IRESD(900000),ISTS(1000000),IENS(1000000),nres,
     1    allins(1000000),nforce,navar,ngvar,ncvar,nuvar,mresd(1000000)
       COMMON /GAMA/ATGC(1000000)
c       common /options/cutang,cuteng,cutoff,hetatm
        common /options/cutang,cuteng,cutoff,hetatm,modelno
       common /residue/resid(900000)
       common /occur/occ(900000)
	distcut=cutoff
	anglcut=cutang
       invdegree=3.14159/180.0 
	rangct=anglcut-20.0
	if(atgc(ists(ires)).eq.ib1) then
          ks=0
          m=0
          do j = ists(ires),iens(ires) 
             if(named(j).eq.at11) ks=j
             if(named(j).eq.atp11)M=j
	  enddo
          if(ks.ne.0.and.m.ne.0)then
          call verify(ires,ks,M,at21,ib2,dist1,thtai1,NN1,atmN,
     1    distcut,0,anglcut,0,0)
          endif
          if(nn1.ne.0) then
            ks=0
            m=0
            do k =ists(ires),iens(ires)
              if(named(k).eq.at12)KS=k
              if(named(k).eq.atp12)M=k
            enddo
          if(ks.ne.0.and.m.ne.0)then
            call verify(ires,KS,M,at22,ib2,dist2,thtai2,NN2,atmN,
     1   distcut,nn1,anglcut,0,0)
          endif
          ks=0
          m=0
              do k = ists(ires),iens(ires)
               if(named(k).eq.at13)KS=k
               if(named(k).eq.atp13)M=k
              enddo
          if(ks.ne.0.and.m.ne.0)then
              call verify(ires,KS,M,at23,ib2,dist3,thtai3,NN3,atmN,
     1   distcut,nn1,anglcut,0,0)
            endif
c  
c  Backtrace loop begins
c
          if(nn2.ne.0.or.nn3.ne.0)then
          ks=0
          m=0
            do k = ists(nn1),iens(nn1)
              if(named(k).eq.at21)KS=k
              if(named(k).eq.atp21)M=k
            enddo
          if(ks.ne.0.and.m.ne.0)then
            call verify(nn1,KS,M,at11,ib1,dist1,thtai4,NN4,atmN,
     1   distcut,ires,rangct,0,0)
            endif
            KS11=KS
            KS12=atmN
          ks=0
          m=0
            do k = ists(nn1),iens(nn1)
               if(named(k).eq.at22)KS=k
               if(named(k).eq.atp22)M=k
            enddo
          if(ks.ne.0.and.m.ne.0)then
            call verify(nn1,KS,M,at12,ib1,dist2,thtai5,NN5,atmN,
     1   distcut,ires,rangct,0,0)
            endif
            KS21=KS
            KS22=atmN
c            if(NN4.eq.0)then
          ks=0
          m=0
              do k = ists(nn1),iens(nn1)
                if(named(k).eq.at23)KS=k
                if(named(k).eq.atp23)M=k
              enddo
c
          if(ks.ne.0.and.m.ne.0)then
              call verify(nn1,KS,M,at13,ib1,dist3,thtai6,NN6,atmN,
     1   distcut,ires,rangct,0,0)
          endif
              KS31=KS
              KS32=atmN
c            endif
            calenrg12=(dist1-3.0)**2+(dist2-3.0)**2 + 0.5*(
     1      ((thtai1-180.0)*invdegree)**2+((thtai2-180.0)*invdegree)
     2       **2+((thtai4-180.0)*invdegree)**2+
     3       ((thtai5-180.0)*invdegree)**2)
c
            calenrg13=(dist1-3.0)**2+(dist3-3.0)**2 + 0.5*(
     1      ((thtai1-180.0)*invdegree)**2+((thtai3-180.0)*invdegree)
     2       **2+((thtai4-180.0)*invdegree)**2+
     3       ((thtai6-180.0)*invdegree)**2)
c
             if(calenrg12.lt.calenrg13)then
              calenrg = calenrg12
              KS21=KS21
              KS22=KS22
             else
              calenrg = calenrg13
              KS21=KS31
              KS22=KS32
             endif
c
c	bhatta for debug April 18, 2006
c	write(*,*) 'Ener',dist1,dist2,thtai1,thtai2,thtai3,thtai4,thtai5,
c     1    thtai6,calenrg12, calenrg13,calenrg
              call orientation(KS11,KS12,KS21,KS22,config,anglconv)
c
	    ncback=nn4 + nn5 + nn6
            if(ncback.ge.2*ires) then
             if(calenrg.le.cuteng.and.anglconv.le.90.0)then
c
               base(1,ires) = base(1,ires)+1
               base(base(1,ires)+1,ires) = nn1
	       feature(base(1,ires),ires) = hbftr
	       type(base(1,ires),ires) = config 
               energy(base(1,ires),ires) = calenrg
              yangle(base(1,ires),ires)=anglconv
c
               base(1,nn1) = base(1,nn1)+1
               base(base(1,nn1)+1,nn1) = ires
	       feature(base(1,nn1),nn1) = hbftrbck
	       type(base(1,nn1),nn1) = config 
               energy(base(1,nn1),nn1) = calenrg
	call orientation(KS21,KS22,KS11,KS12,config,anglconv)
               yangle(base(1,nn1),nn1)=anglconv
            endif
            endif
          endif 
          endif 
	endif
        return
        end
C
C	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
       subroutine findpair(ires,ib1,atp11,at11,at21,atp21,ib2,
     1   atp12,at12,at22,atp22,config,hbftr,hbftrbck,nsg)
       EXTERNAL VERIFY
       integer ib1,ib2,ATGC,KS,NN1,NN2,ires,atmN,base,NN3,NN4,KS11,KS12
     1          ,KS21,KS22
       real   distcut,anglcut,energy,calenrg,anglconv,yangle,occ
	real invdegree
       character*4 at11,atp11,at21,at12,atp12,at22,named,atp21,
     1              atp22
       character*3 hbftr,feature,hbftrbck,resid,resd
       character*1 config,type,allins
       character*4 chaind
       COMMON /HBQUA/BASE(21,20000),energy(20,20000),
     1                    yangle(20,20000)
       COMMON /HBQTYP/TYPE(20,20000),FEATURE(20,20000)
c       COMMON /NAMES/NAMED(900000),RESD(900000)
       COMMON /NAMES/NAMED(900000),RESD(900000)
       COMMON /CHAINS/CHAIND(900000)
       COMMON /ALPHA/XD(900000),YD(900000),ZD(900000)
       COMMON /BETA/IRESD(900000),ISTS(1000000),IENS(1000000),nres,
     1     allins(1000000),nforce,navar,ngvar,ncvar,nuvar,mresd(1000000)
       COMMON /GAMA/ATGC(1000000)
c       common /options/cutang,cuteng,cutoff,hetatm
        common /options/cutang,cuteng,cutoff,hetatm,modelno
       common /residue/resid(900000)
       common /occur/occ(900000)
	ksugar=0
	if(nsg.eq.0) then
	distcut=cutoff
	else
	distcut=cutoff+2.0
	endif
	anglcut=cutang
       if(hbftr(1:1).eq.'+'.or.hbftrbck(1:1).eq.'+'.or.hbftr(1:1).eq.'z'
     1 .or.hbftrbck(1:1).eq.'z')anglcut=anglcut+30.0
C  More strict evaluation of PROTONATED base pair
C  bhatta, Jan 2010
	angcut=anglcut - 20.0
        invdegree=3.14159/180.0 
	if(ATGC(ists(ires)).eq.ib1)then
         ks=0
          m=0
          do j = ists(ires),iens(ires)  
           if(named(j).eq.at11) ks=j
           if(named(j).eq.atp11)M=j
          enddo
          if(ks.ne.0.and.m.ne.0)then
	  call verify(ires,KS,M,at21,ib2,dist1,thtai1,NN1,atmN,
     1     cutoff,0,anglcut,0,0)
          endif
	  if(nn1.ne.0) then
         ks=0
          m=0
               do k = ists(ires),iens(ires)
	         if(named(k).eq.at12)KS = k
                 if(named(k).eq.atp12)M = k
               enddo
	  if(nsg.ne.0) then
	    do k=ists(ires),iens(ires)
	      if(named(k).eq.'O2* ') ksugar=k
	    enddo
	  endif
c
c   In this case distance calculation is done between ksugar and at22
c   while linearity check is done with at12, atp12 and at22
c
          if(ks.ne.0.and.m.ne.0)then
               call verify(ires,KS,M,at22,ib2,dist2,thtai2,NN2,atmN,
     1       distcut,nn1,anglcut,ksugar,0)
          endif
             if(nn1.eq.nn2)then
         ks=0
          m=0
c               anglcut=anglcut-20.0
               
c
               do k = ists(nn1),iens(nn1)
	         if(named(k).eq.at21)KS = k
                 if(named(k).eq.atp21)M=k
               enddo
          if(ks.ne.0.and.m.ne.0)then
	       call verify(nn1,KS,M,at11,ib1,dist3,thtai3,NN3,atmN,
     1    cutoff,ires,angcut,0,0)
          endif
             KS11=KS
             KS12=atmN
         ks=0
          m=0
               do k = ists(nn1),iens(nn1)
	           if(named(k).eq.at22)KS = k
                   if(named(k).eq.atp22)M = k
               enddo
          if(ks.ne.0.and.m.ne.0)then
	    if(nsg.eq.0) then
               call verify(nn1,KS,M,at12,ib1,dist4,thtai4,NN4,atmN,
     1      distcut,ires,angcut,0,0)
	    else
c  
c   In this case distance calculation is done betn. at22 & O2* while
c   linearity check is done using at22, atp22 and at12
c
               call verify(nn1,KS,M,at12,ib1,dist4,thtai4,NN4,atmN,
     1      distcut,ires,angcut,0,1)
            endif
          endif
             KS21=KS
             KS22=atmN
c 
             calenrg=(dist1-3.0)**2+(dist2-3.0)**2 + 0.5*(
     1      ((thtai1-180.0)*invdegree)**2+((thtai2-180.0)*invdegree)
     2       **2+((thtai3-180.0)*invdegree)**2+
     3       ((thtai4-180.0)*invdegree)**2)

             call orientation(KS11,KS12,KS21,KS22,config,anglconv)

               if((nn3.eq.nn4).and.(nn3.eq.ires))then
                if(calenrg.le.cuteng.and.anglconv.le.90.0)then
c               
                  base(1,ires) = base(1,ires)+1
                  base(base(1,ires)+1,ires) = nn1
	          feature(base(1,ires),ires) = hbftr
	          type(base(1,ires),ires) = config
       	          energy(base(1,ires),ires)=calenrg
                  yangle(base(1,ires),ires)=anglconv
c
                  base(1,nn1) = base(1,nn1)+1
                  base(base(1,nn1)+1,nn1) = ires 
	          feature(base(1,nn1),nn1) = hbftrbck
	          type(base(1,nn1),nn1) = config
                  energy(base(1,nn1),nn1)=calenrg
c                  call orientation(KS21,KS22,KS11,KS12,config,anglconv)
                  yangle(base(1,nn1),nn1)=anglconv
                endif
               endif  
	     endif
          endif
	endif
        return
        end

c       ###############################################################

        subroutine verify(srcres,srcatm,atmprcsr,trgtatm,base2,
     1    optdist,thetai,NN,atmN,cut,nfind,thetamin,nsug,ksug)
c
        integer base2,ATGC,atmN,srcres,srcatm,atmprcsr,nfind,NN
	real optdist,cut,thetai,thetamin,occ       
        CHARACTER*4 trgtatm,named,chaind
        character*3 resid,resd
        character*1 allins
        COMMON /ALPHA/XD(900000),YD(900000),ZD(900000)
	COMMON /BETA/IRESD(900000),ISTS(1000000),IENS(1000000),nres,
     1    allins(1000000),nforce,navar,ngvar,ncvar,nuvar,mresd(1000000)
        COMMON /GAMA/ATGC(1000000) 
        COMMON /NAMES/NAMED(900000),RESD(900000)
        COMMON /CHAINS/CHAIND(900000)
        common /residue/resid(900000)
        common /occur/occ(900000)
	degree = 180.0/3.1459
        NN = 0
	distsg = 0.0
        atmN = 0
        optdist = 1000.0
        thetai = 0.0
	if(nfind.eq.0) then
	   kbegn=1
	   kend=nres
	else
	   kbegn=nfind
	   kend=nfind
	endif
	do k = kbegn,kend
         if(ATGC(ists(k)).eq.base2)then
         if((chaind(ists(k)).ne.chaind(ists(srcres))).or.
     1      (abs(k-srcres).gt.2).or.(nforce.eq.1.and.k.ne.srcres))then
	   if(ksug.eq.1) then
	    do l=ists(k),iens(k)
	      if(named(l).eq.'O2* ') then
              distsg=(xd(srcatm)-xd(l))**2+(yd(srcatm)-yd(l))**2+
     1              (zd(srcatm)-zd(l))**2
	      endif
	    enddo
	   endif
	      
            do l = ists(k),iens(k)
	     if(named(l).eq.trgtatm)then
	      if(nsug.eq.0) then
              distance=(xd(srcatm)-xd(l))**2+(yd(srcatm)-yd(l))**2+
     1              (zd(srcatm)-zd(l))**2
	     else
              distance=(xd(nsug)-xd(l))**2+(yd(nsug)-yd(l))**2+
     1              (zd(nsug)-zd(l))**2
	     endif
	      
              if(distance.le.cut.or.(ksug.eq.1.and.distsg.le.cut))then
                bonddist=sqrt(distance)
                 call linearity(srcatm,atmprcsr,l,theta)
                  thetad = theta*degree
                  if(thetad.ge.thetamin)then
                    if(NN.ne.0)then
                      if(thetad.gt.thetai)then
                         NN = k
                         atmN = l
                         optdist = bonddist
                         thetai = thetad
                      endif   
                    else 
                         NN = k
                         atmN = l
                         optdist = bonddist
                         thetai = thetad
	            endif
	          endif
              endif  
             endif
            enddo
          endif
         endif
        enddo
c	bhatta for debugging April 18, 2006
c	write(*,*) cut,optdist,thetamin,thetai
        return
        end


c       -----------------------------------------------------------#####

c        check for linearity of hydrogen bonds betn. a basepair
c       ----------------------------------------------------------####

        subroutine linearity(j,M,l,theta)
        COMMON /ALPHA/XD(900000),YD(900000),ZD(900000)

            x1=(xd(j)-xd(M))
            y1=(yd(j)-yd(M))
            z1=(zd(j)-zd(M))

            x2=(xd(j)-xd(l))
            y2=(yd(j)-yd(l))
            z2=(zd(j)-zd(l))

	    anum = x2*x1 + y2*y1 + z2*z1
	    dist1 = sqrt(x1*x1 + y1*y1 + z1*z1)
	    dist2 = sqrt(x2*x2 + y2*y2 + z2*z2)
	    theta = acos (anum/(dist1*dist2))
c        theta = Acos(((x2*x1)+(y2*y1)+(z2*z1))/((sqrt(x1**2.0+
c     1          y1**2.0+z1**2.0))*(sqrt(x2**2.0+y2**2.0+z2**2.0))))


        return
        end


c       #################################################################


c      -----------------------------------------------------------######
c      Assigning the orientation of the basepairs [ Cis or Trans ]
c      ----------------------------------------------------------#######

	subroutine orientation(atm11,atm12,atm21,atm22,orient,anglconv)
        real angle,anglconv
        integer atm11,atm12,atm21,atm22
        character*1 orient
        COMMON /ALPHA/XD(900000),YD(900000),ZD(900000)
c
        degree = 180.0/3.1459
	if(atm11.ne.0.and.atm12.ne.0.and.atm21.ne.0.and.atm22.ne.0) then
		! bhatta March 17, 2021
c	write(6,*) 'atm11, atm21, atm12, atm22',atm11,atm21,atm12,atm22
c
        x11 = xd(atm21)-xd(atm11)
        y11 = yd(atm21)-yd(atm11)
        z11 = zd(atm21)-zd(atm11)
c
        x22 = xd(atm22)-xd(atm12)
        y22 = yd(atm22)-yd(atm12)
        z22 = zd(atm22)-zd(atm12)
c        if(orient.eq.'T')then
c         x11 = -x11
c         y11 = -y11
c         z11 = -z11
c        
c         x22 = x22
c         y22 = y22
c         z22 = z22
c        else
c         x11 = x11
c         y11 = y11
c         z11 = z11
c        
c         x22 = x22
c         y22 = y22
c         z22 = z22
c       endif
         

        angle = Acos(((x22*x11)+(y22*y11)+(z22*z11))/((sqrt(x11**2.0+
     1       y11**2.0+z11**2.0))*(sqrt(x22**2.0+y22**2.0+z22**2.0))))
       
        anglconv = angle*degree
c        if(anglconv.ge.0.and.anglconv.le.30.0)then
c         orient = 'C'
c        elseif(anglconv.ge.150.0.and.anglconv.le.180.0)then
c         orient = 'T'
c        endif
	endif
       return
       end
C     ************** To clean PDBFILE ****************************
C
	subroutine cleanpdb(filenm,nocor,nofasta)
       character*3 resid,resd
       CHARACTER*4 NAMED,atomnm(9)
       character*1 fastaseq(100000),allins
       CHARACTER*4 CHAIND,pcd,pchaind,molid
       integer presd, prd,nstfrg
        COMMON /NAMES/NAMED(900000),RESD(900000)
       COMMON /CHAINS/CHAIND(900000)
       COMMON /ALPHA/XD(900000),YD(900000),ZD(900000)
       COMMON /BETA/IRESD(900000),ISTS(1000000),IENS(1000000),nres,
     1     allins(1000000),nforce,navar,ngvar,ncvar,nuvar,mresd(1000000)
       common /residue/resid(900000)
       common /occur/occ(900000)
       common /num/kresd(1000000),prd(1000000)
       common /cum/pchaind(1000000),pcd(1000000),molid

       real occ
       CHARACTER*80 filenm
	data atomnm/'N1  ','C2  ','N3  ','C4  ','C5  ','C6  ','N7  ',
     1 'C8  ','N9  '/
	
	numresb=nres
        ild=nres
834    format(i5,i5,a2)
        do i=1,nres
         prd(i)=kresd(i)
         pcd(i)=pchaind(i)
c         write(6,*)i,prd(i),kresd(i),pcd(i),pchaind(i)
        end do
	i=1
	do while (i.le.numresb)
	  natominr=0
	  ncafnd=0
	  do k=ists(i),iens(i)
	    if(named(k)(1:1).eq.' ')then
	      named(k)(1:3)=named(k)(2:4)
	      named(k)(4:4)=' '
	    endif
	    if(named(k)(1:1).eq.' ')then
	      named(k)(1:3)=named(k)(2:4)
	      named(k)(4:4)=' '
	    endif
            if(named(k)(3:3).eq.'''')then
               named(k)(3:3)='*'
            endif
	    do m=1,9
	      if(atomnm(m).eq.named(k)) natominr=natominr+1
	    enddo
	    if(named(k).eq.'C1* ') ncafnd=1
	  enddo
c        write(6,*) 'MolID=',molid
c	  if(molid.ne.' '.and.molid.ne.pcd(i)) then
c            write(6,*) 'NatomInR',natominr,'nCAfnd',ncafnd,i,numresb
	  if(natominr.lt.6.or.ncafnd.eq.0) then
	    do m=i,nres-1
	      ists(m)=ists(m+1)
	      iens(m)=iens(m+1)
              prd(m)=prd(m+1)
              pcd(m)=pcd(m+1)
	    enddo
	    numresb=numresb-1
	    i=i-1
	  endif
c   	  endif
	  i=i+1
	enddo 
c	write(52,75)nres
	write(52,76)numresb
	nres=numresb
75	format('#HEADER   Original number of residues in structure ',
     1  'file: ',i6)
76	format('#HEADER   Cleaned number of residues: ',i6)

        nstfrg=1
        do i=1,nres
           fastaseq(i)=resd(ists(i))(3:3)
	   if(ists(i).gt.1) then	! bhatta March 17, 2021
           if((i.gt.1.and.(chaind(ists(i)).ne.chaind(ists(i)-1))).or.i
     1     .eq.nres) then
              nrsfrg=i-1
              if(i.eq.nres) nrsfrg=i
	      if(nofasta.eq.0) then
                write(77,114) filenm(1:20),chaind(ists(i)-1)
                write(77,113) (fastaseq(j),j=nstfrg,nrsfrg)
	      endif
              nstfrg=i
           endif
	   endif	! bhatta March 17, 2021
           do l=ists(i),iens(i)
             if(ists(i).gt.1) then 	! bhatta March 17, 2021
              if(named(l).eq.named(l-1).and.resd(l).eq.resd(l-1).and.
     1  xd(l).eq.xd(l-1).and.yd(l).eq.yd(l-1).and.zd(l).eq.zd(l-1)) then
	        if(nocor.eq.0) then
                  write(9,112)
	        endif
                return
              endif
             endif	! bhatta March 17, 2021
	      if(nocor.eq.0) then
              write(9,111)l,named(l),resd(l),i,xd(l),yd(l),zd(l),occ(l)
	      endif
c              prd(l)=i
C Dirty Fix
          end do
        end do
        do i=1,nres
         kresd(i) = prd(i)
         pchaind(i) = pcd(i)
c         write(6,*)i,prd(i),kresd(i),pcd(i),pchaind(i)
        end do

	if(nocor.eq.0) then
	  write(9,112)
	endif
111     format('ATOM',1X,I6,2X,A4,A3,I6,4X,3F8.3,2x,f4.2,F6.2)
112	format('END')
113     format(80a1)
114     format('>',a20,':',a4)
	close(unit=9)
	close(unit=77)
	return
	end
	subroutine showoptn(ierror)
	write(6,1)
	write(6,2)
	write(6,3)
	write(6,4)
	write(6,10)
	write(6,5)
	write(6,6)
	write(6,7)
	write(6,8)
	write(6,9)
	write(6,11)
        write(6,13)
	write(6,12)
1     format('The following options can be used:')
2     format('-HD [value] to set default hydrogen bond distance ',
     1 'cutoff (default = 3.8)')
3     format('-VA [value] to set default pseudo angle cutoff ',
     1 '(default = 120.0)')
4     format('-EN [value] to set default E-value cutoff ',
     1 '(default = 1.8)')
5     format('-HT to include HETATM entries in PDB')
6     format('-CH to avoid identification of base pairs stabilized', 
     1 ' by C-H...O/N H-bonds') 
7     format('-SG to avoid identification of base pairs involving', 
     1 'sugar O2'' atoms')
8     format('-AB to avoid base pairing between residue no. i and',
     1 ' i+1') 
9     format('-OL to avoid printing base pairing information w.r.t.',
     1  'the second strand.  This is suitable for simple ',
     1  'oligonucleotides')
10    format('-ML [character] to select desired chain identifier ',
     1  'in PDB File (default = all)')
11    format('-NMR to calculate base pairing information for the first',
     1  ' model NMR derived structure')
13    format('-MD [number] to calculate base pairing information for',
     1  'a particular NMR model [number]')
12    format('-CIF to calculate base pairing information from mmCIF'
     1  ' formatted file (PDB is the default). This considers only the',
     1  ' first model of NMR derived structures')
	ierror=1
	return
	end


C TO FIND OUT THE LENGTH OF THE LARGEST CONTINUOUS HELIX IN RNA 
C THE PROGRAM USES *.out FILES AS INPUT
C AUG 8; 4:10

C LRNA : TOTAL LENGTH OF RNA
C IHLX(X,Y) : INFORMATION OF ANY HELIX AT LEAST 3 BP IN LENGTH
C             - HELIX STARTING POINT AND LENGTH
C NHLX : NO. OF HELICES GREATER THAN 3 BP             
C NWC : NO.OF HELICES CONTAINING ALL WATSON-CRICK BP
C NNW : NO.OF HELICES CONTAINING AT LEST ONE NON-WATSON-CRICK BP
       subroutine rnahelix(nmpass,filenm,nnf,modelno,nocsv,nodat,nodbn,
     1    nohlx)
       DIMENSION IRES1(1000000),IRES2(1000000),
     1       nrstk(999,2),ndoub(999,2)
       DIMENSION  npdbres(1000000)
       dimension strdbn(1000000)
       INTEGER nchnb(200),nchne(200)
       CHARACTER*1 BASE1(1000000),BASE2(1000000),strseq,strdbn,
     1          allins
       CHARACTER*80 nmpass,filenm,file2,text2wr
      	character*132 line
       CHARACTER*5 BPAIR(1000000)
       character*15 prnvar,prnvar1
       character*3 resd,bs1,bs2
       character*4 pcd,pchaind,namechn(500),named,molid
       common /str/strseq(1000000)
       common /cum/pchaind(1000000),pcd(1000000),molid
        COMMON /NAMES/NAMED(900000),RESD(900000)
       COMMON /BETA/IRESD(900000),ISTS(1000000),IENS(1000000),nres,
     1     allins(1000000),nforce,navar,ngvar,ncvar,nuvar,mresd(1000000)

       OPEN(UNIT=1,FILE=nmpass,STATUS='OLD',ERR=100)
c       OPEN(UNIT=2,FILE='nwhelix.dat')
c       OPEN(UNIT=3,FILE='wchelix.dat')
c       open(unit=8,file='pshelix.dat')
c	write(6,*) 'In RNAHELIX subroutine'
       nn=index(nmpass,'.')
        file2=nmpass
c	if(nohlx.eq.0) then
          file2(nn:nn+3)='.hlx'
          open(unit=9,file=file2)
c	endif
       
        text2wr=nmpass(1:index(nmpass,'.')-1)
         lpstart=index(text2wr,'/')
       do while(lpstart.gt.0)
         lpstart=index(text2wr,'/')
         if(lpstart.ne.0)  text2wr=text2wr(lpstart+1:80)
       enddo 
       
       NNW=0
       NWC=0
       NStkHlx=0
       k=0
       do while (k.eq.0)
         READ(1,5,END=20)LINE
         IF(LINE(65:74).EQ.'Equivalent') THEN
           DO I=1,1000000
             read(1,14,END=20) line
           enddo
         endif
       enddo
20     continue
       lrna=i-1
C *********************************************************************       
       close(unit=1)
       OPEN(UNIT=1,FILE=nmpass,STATUS='OLD',ERR=100)
           write(81,170) text2wr(1:index(text2wr,' ')),lrna
c        write(81,164) nmpass(1:4),lrna
       K=0
       DO WHILE(K.EQ.0)
         READ(1,5,END=10)LINE
         IF(LINE(65:74).EQ.'Equivalent') THEN
           DO I=1,1000000
             read(1,14,END=10) line
c           READ(line,15)IRES1(I),npdbres(i),BASE1(I),IRES2(I),
c     1   BASE2(I),BPAIR(I)
           READ(line,15)IRES1(I),npdbres(i),BS1,IRES2(I),
     1   BS2,BPAIR(I),ires3,ires4
14    format(a132)
15    FORMAT(I6,i8,t16,A3,7X,I6,t39,A3,t50,A5,t63,i6,t100,i6)
        call renbase(bs1,base1(i))
        call renbase(bs2,base2(i))
c           write(6,15) ires1(i),(i),base1(i),ires2(i),base2(i),
c     1 bpair(i)
           write(81,162) ires1(i),base1(i),ires2(i),ires3,ires4
           ENDDO
         K=1
         ENDIF
       ENDDO

10     CONTINUE
       LRNA=I-1
	do kres=1,lrna
         strdbn(kres)='.'
        enddo
	do kres=1,lrna
         if(ires2(kres).ne.0) then
          if(kres.eq.(ires2(ires2(kres))).and.ires2(kres).gt.kres)then
c           if(BPAIR(kres)(1:4).eq.'W:WC'.and.strdbn(kres).ne.')')then
           if(strdbn(kres).ne.')')then
c            if(ires2(kres-1).eq.ires2(kres)+1.and.ires2(kres+1).eq.
c     1  ires2(kres)-1) then
C
C  Searching for Pseudo-knot
C
             noopen=0
             noclose=0
             do kc=kres,ires2(kres)
               if(strdbn(kc).eq.'(') noopen=noopen+1
               if(strdbn(kc).eq.')') noclose=noclose+1
             enddo
c        write(6,*) 'kres=',kres,noopen,noclose
             if(noopen.lt.noclose) then
               noop2=0
               nocl2=0
               do kc=kres,ires2(kres)
                if(strdbn(kc).eq.'[') noop2=noop2+1
                if(strdbn(kc).eq.']') nocl2=nocl2+1
               enddo
               if(noop2.lt.nocl2) then
               strdbn(kres)='{'
               strdbn(ires2(kres))='}'
               else
               strdbn(kres)='['
               strdbn(ires2(kres))=']'
               endif
             else 
               strdbn(kres)='('
               strdbn(ires2(kres))=')'
             endif
c        write(6,*)'WC',kres,ires2(kres),strdbn(kres),strdbn(ires2(kres))
c            endif
           else
c             strdbn(kres)='('
c             strdbn(ires2(kres))=')'
c        write(6,*) 'Found Non-WC pair',kres,ires2(kres)
           endif
          endif
         endif
        enddo
c	  write(23,168) (strseq(i),i=1,nres)

C *********************************************************************
C Finding Pseudo Continuous Helix 
C
       k=0
       DO J=1,(lrna-6)
        IF((IRES2(J).NE.0).and.
     5    (IRES2(J+3).NE.(IRES2(J+2)-1).AND.(IRES2(J+3).NE.0)).and.
     1       (IRES2(J+1).EQ.(IRES2(J)-1).AND.(IRES2(J+1).NE.0)).and. 
     2       (IRES2(J+2).EQ.(IRES2(J+1)-1).AND.(IRES2(J+2).NE.0)).and.
     3       (IRES2(J+4).EQ.(IRES2(J+3)-1).AND.(IRES2(J+4).NE.0)).and.
     4       (IRES2(J+5).EQ.(IRES2(J+4)-1).AND.(IRES2(J+5).NE.0))) THEN
           m=0
           lext=j-1
           NStkHlx=NStkHlx+1
           nrstk(NStkHlx,1)=IRES1(J)
           do while(ires2(lext).ne.0.and.(ires2(lext).eq.
     1                 ires2(lext+1)+1))
C
c Extension of minimal Pseudohelix towards 5'-end side of first strand'
C
             nrstk(NStkHlx,1)=IRES1(lext)
             lext=lext-1
             enddo
	     lext=j+6
	do while(ires2(lext).ne.0.and.(ires2(lext).eq.ires2(lext-1)-1))
C
c Extension of minimal Pseudohelix towards 3'-end side of first strand'
C
               lext=lext+1
	     enddo
	  nrstk(NstkHlx,2)=IRES1(lext-1)
        ENDIF
       ENDDO
       kcontin=k
c	write(*,*) 'NStkHlx',NstkHlx
C
C Now finding out the normal antiparallel double helical regions
C which do not appear within Stacked PseudoContinuous Helices
C
       j=1
       ndef=0
       NdbHlx=0
c	write(*,*) 'StkHlx',(nrstk(jl,1),nrstk(jl,2),jl=1,NstkHlx)
       do while(j.le.(lrna-3))
         if(Ires2(j).gt.Ires1(j)) then
           if((Ires2(j+1).eq.(Ires2(j)-1)).and.(ires2(j+2).eq.
     1                                       (ires2(j+1)-1))) then
             NdbHlx=NdbHlx+1
             ndoub(NdbHlx,1)=ires1(j)
             j=j+2
             ndoub(NdbHlx,2)=ires1(j)
             do while(j.le.lrna)
               if(Ires2(j+1).eq.0) exit
               if(Ires2(j+1).eq.(Ires2(j)-1)) then
                 ndoub(NdbHlx,2)=ires1(j+1)
                 j=j+1 
	       else
c	       write(*,*) 'Exiting DblHlxFinding',NdbHlx,j, ires2(j)
	         exit
               endif
             enddo
c	       write(*,*) 'Exited DblHlxFinding',NdbHlx,j, ires2(j)
c
c   Finding Hairpin loops
c
           if(ires2(j).eq.ires1(j+3).and.ires2(j+1).eq.0.and.ires2(j+2)
     1         .eq.0) then
c            write(*,*) 'Diloop found at',ires1(j),ires2(j)
            strseq(ires1(j+1))='L'
            strseq(ires1(j+2))='L'
           elseif(ires2(j).eq.ires1(j+4).and.ires2(j+1).eq.0.and.
     1         ires2(j+2).eq.0.and.ires2(j+3).eq.0) then
c            write(*,*) 'Triloop found at',ires1(j),ires2(j)
            strseq(ires1(j+1))='L'
            strseq(ires1(j+2))='L'
            strseq(ires1(j+3))='L'
           elseif(ires2(j).eq.ires1(j+5).and.ires2(j+1).eq.0.and.
     1     ires2(j+2).eq.0.and.ires2(j+3).eq.0.and.ires2(j+4).eq.0) then
c            write(*,*) 'Tetraloop found at',ires1(j),ires2(j)
            strseq(ires1(j+1))='L'
            strseq(ires1(j+2))='L'
            strseq(ires1(j+3))='L'
            strseq(ires1(j+4))='L'
           elseif(ires2(j).eq.ires1(j+6).and.ires2(j+1).eq.0.and.
     1     ires2(j+2).eq.0.and.ires2(j+3).eq.0.and.ires2(j+4).eq.0.
     1     and.ires2(j+5).eq.0) then
c            write(*,*) 'Pentaloop found at', ires1(j),ires2(j)
            strseq(ires1(j+1))='L'
            strseq(ires1(j+2))='L'
            strseq(ires1(j+3))='L'
            strseq(ires1(j+4))='L'
            strseq(ires1(j+5))='L'
           elseif(ires2(j).eq.ires1(j+7).and.ires2(j+1).eq.0.and.
     1     ires2(j+2).eq.0.and.ires2(j+3).eq.0.and.ires2(j+4).eq.0.
     1     and.ires2(j+5).eq.0.and.ires2(j+6).eq.0) then
            strseq(ires1(j+1))='L'
            strseq(ires1(j+2))='L'
            strseq(ires1(j+3))='L'
            strseq(ires1(j+4))='L'
            strseq(ires1(j+5))='L'
            strseq(ires1(j+6))='L'
           elseif(ires2(j).eq.ires1(j+8).and.ires2(j+1).eq.0.and.
     1     ires2(j+2).eq.0.and.ires2(j+3).eq.0.and.ires2(j+4).eq.0.
     1     and.ires2(j+5).eq.0.and.ires2(j+6).eq.0.and.ires2(j+7).eq.0
     1     ) then
            strseq(ires1(j+1))='L'
            strseq(ires1(j+2))='L'
            strseq(ires1(j+3))='L'
            strseq(ires1(j+4))='L'
            strseq(ires1(j+5))='L'
            strseq(ires1(j+6))='L'
            strseq(ires1(j+7))='L'
           endif
C Considers upto 7 residue long hairpin loops now. bhatta March 2020

           endif
            
c	 enddo
C
C 
C   Checking whether the same helix was earlier detected as part of a continuous helix
C
c           ndef=0
c	if(j.eq.ndoub(ndbhlx,2)) then
c	write(*,*) 'Going to find detected continuous helix',j

!         do k=1,NstkHlx
!           if((j.le.nrstk(k,2).and.j.ge.nrstk(k,1)).or.(Ires2(j).le.
!     1         nrstk(k,2).and.Ires2(j).ge.nrstk(k,1))) then
!	     write(*,*) j,' already registrared Helix',NdbHlx,j,Ires2(j),nrstk(k,1)
!             j=nrstk(k,2)+1
c	     ndef=1
!             NdbHlx=NdbHlx-1
!             exit
!           endif
!         enddo
c	 endif
         endif
         j=j+1
         enddo
c       write(*,*) 'NdbHlx',(ndoub(jl,1),ndoub(jl,2),jl=1,NdbHlx)

       do kl=1,NstkHlx
         do j=nrstk(kl,1),nrstk(kl,2)
	   if(nohlx.eq.0) then
           write(9,17) ires1(j),ires2(j),bpair(j)(1:3),bpair(j)(4:4),
     1  base1(j),base2(j)
	   endif
           if(strseq(ires1(j)).ne.'N'.and.strseq(ires2(j)).ne.'N')then
             if(strseq(ires1(j)).ne.'T') strseq(ires1(j))='H'
             if(strseq(ires2(j)).ne.'T') strseq(ires2(j))='H'
	   endif
         enddo
	 if(nohlx.eq.0) then
         write(9,17) 0
	 endif
       enddo
	if(nohlx.eq.0) then
       write(9,17) 0
	endif
       do kl=1,NdbHlx
         do j=ndoub(kl,1),ndoub(kl,2)
	   if(nohlx.eq.0) then
           write(9,17) ires1(j),ires2(j),bpair(j)(1:3),bpair(j)(4:4),
     1   base1(j),base2(j)
	   endif
           if(strseq(ires1(j)).ne.'N'.and.strseq(ires2(j)).ne.'N') then
C Added the conditions so that Triplets are not considered Helix. bhatta
C May 9, 2018
             if(strseq(ires1(j)).ne.'T') strseq(ires1(j))='H'
             if(strseq(ires2(j)).ne.'T') strseq(ires2(j))='H'
           endif
         enddo
        if(ndoub(kl+1,1).gt.0.and.(ndoub(kl+1,1)-ndoub(kl,2)).le.5) then
c	write(6,*) ndoub(kl+1,1),ndoub(kl,2)
           if(ires2(ndoub(kl+1,1))-ires2(ndoub(kl,2)).eq.-1) then
c             write(*,*) 'possible bulge at', ndoub(kl,2) 
             do kk=ndoub(kl,2)+1,ndoub(kl+1,1)-1
              strseq(kk)='B'
             enddo
           endif
         endif
	 if(nohlx.eq.0) then
         write(9,17) 0
         endif
       enddo

c	write(*,*) 'M1 - size of helices',(m1(l),l=1,nhlx)
c	write(*,*) 'IHLX(1),IHLX(2)',(ihlx(l,1),ihlx(l,2),l=1,nhlx)

c      WRITE(*,*)'Total number of helices:',(NNW+NWC)
c      WRITE(*,*)'No. of antiparallel helices with all W-C base pairs:
c     1 ',NWC
c      WRITE(*,*)'No. of antiparallel helices with at least one non-W-C
c     1 base pairs:',NNW
      close(unit=1) 
       nchain=1
       nchnb(nchain)=1
       namechn(nchain)=pchaind(1)
c      write(*,*) (pchaind(kk),kk=1,nres)
      do i=2,nres
        if(pchaind(i).ne.pchaind(i-1)) then
c           write(*,*)'Chain IDs change at',i
           nchne(nchain)=i-1
           if(nchain.le.499) then
              nchain=nchain+1
           else
              write(6,*) 'There are more than 100 chains in the pdb'
	      ierror=1
              return
           endif
           nchnb(nchain)=i
           namechn(nchain)=pchaind(i)
        endif 
      enddo
      nchne(nchain)=nres
      do ii=1,nchain
c        iff=index(namechn(ii),' ')
c        namechn(ii)=namechn(ii)(1:iff)
        krun=0
	if(nocsv.eq.0) then
        do nche=nchnb(ii),nchne(ii)
          write(prnvar1,*) ires1(nche)
          prnvar1=adjustl(prnvar1)
          iff1=index(prnvar1,' ')-1
          write(prnvar,*)npdbres(nche)
          prnvar=adjustl(prnvar)
          iff=index(prnvar,' ')-1
          iff2=index(namechn(ii),' ')-1
          write(78,*) nmpass(1:4),',',prnvar1(1:iff1),',',prnvar(1:iff)
     1  ,',',base1(nche),',',
     1  allins(ists(nche)),',',namechn(ii)(1:iff2),',',strseq(nche)
        enddo
	endif
c         write(79,169) filenm(1:nnf-1),namechn(ii)
        write(prnvar1,*) modelno
        prnvar1=adjustl(prnvar1)  
        prnvar=adjustl(namechn(ii))
        prnvar=trim(prnvar)
	if(nodat.eq.0) then
         write(79,169) text2wr(1:index(text2wr,' ')),namechn(ii),prnvar1
C         write(6,*) ii,nchnb(ii),nchne(ii), namechn(ii)
         write(79,168) (strseq(i),i=nchnb(ii),nchne(ii))
	endif
c        if(nche.eq.1) then
c         write(80,165) (base1(i),i=nchnb(ii),nchne(ii))
c        else
c          krun=krun+1
c          strdbn(nchne(ii)+1)='&'
c        endif
        
      enddo
      
c        write(80,164) nmpass(1:4),lrna
           write(80,164) text2wr(1:index(text2wr,' ')),lrna
        jwr=1
        do kwr=1,nchne(1)
          strseq(kwr)=base1(kwr)
          jwr=jwr+1
        enddo
        do ii=2,nchain
          strseq(jwr)='&'
          jwr=jwr+1
          do kwr=nchnb(ii),nchne(ii)
            strseq(jwr)=base1(kwr)
            jwr=jwr+1
          enddo
        enddo
          
           write(80,165) (strseq(i),i=1,jwr-1)
        jwr=1
        do kwr=1,nchne(1)
          strseq(kwr)=strdbn(kwr)
          jwr=jwr+1
        enddo
        do ii=2,nchain
          strseq(jwr)='&'
          jwr=jwr+1
          do kwr=nchnb(ii),nchne(ii)
            strseq(jwr)=strdbn(kwr)
            jwr=jwr+1
          enddo
        enddo
          
           write(80,165) (strseq(i),i=1,jwr-1)
5     FORMAT(A80)
C15    FORMAT(I6,7X,A3,5X,I6,7X,A3,6X,A5)
17	format(i5,1x,i5,1x,a3,1x,a1,1x,2a1)
25    FORMAT(I5,1X,I5,A4,1X,A1,1X,'BP',5x,a4)
35    FORMAT(A25)
166   format(a4,1x,i5,a1,1x,a4,1x,a1)
167   format('rm -f ',2a30)
168   format(70a1)
169   format('>',a,':',a,a,' Created by BPFIND')
164   format('>',a,' nts=',i5,' Created by BPFIND')
170   format('#',a,' nts=',i5,' Created by BPFIND')
165   format(20000a1)
163   format('&',20000a1)
162   format(i6,1x,a1,3i6)
100   CONTINUE
        close(unit=1)
        write(line,167)nmpass(1:30),file2(1:30)
c        call system(line)
	close(unit=9)
	close(unit=78)
	close(unit=79)
	close(unit=80)
	close(unit=81)
      return
      END
    
        subroutine readcif (npoorat,modelfound)
        character*1 posei,inscode
        character*4 molid,pcd,pchaind
        character*132 line,type(100)
        character*132 data(5000000)
        character*15 cda(30),atmnm,res,sequ,chain,alternate
        common /cum/pchaind(1000000),pcd(1000000),molid
        real xcoor,ycoor,zcoor,occp
        integer mresd(900000),static
        
        imolid=index(molid,' ')
C        write(6,*) 'MOLID ',molid,imolid
        modelfound=0
        jtype=0
        i=0
        k=0
        numatm=0
        nmrmodel=0
        natmhtm=0
        posei=' '
        inscode=' '
        do while(i.eq.0)
          read(4,1,END=99) line
          if(line(1:5).eq.'loop_') then
            read(4,1) line
            if(line(1:11).eq.'_atom_site.') then
              jtype=jtype+1
              do while(k.eq.0)
!--------------------------
!       if(line(1:20).eq.'_atom_site.group_PDB') natmhtm=jtype
!       if(line(1:24).eq.'_atom_site.label_atom_id') natm=jtype
!       if(line(1:24).eq.'_atom_site.label_comp_id') inres=jtype
!C       if(line)(1:26).eq.'_atom_site.label_entity_id') nseq=jtype
!C       if(line(1:24).eq.'_atom_site.label_asym_id') nchn=jtype
!       if(line(1:23).eq.'_atom_site.auth_asym_id') nchn=jtype
!       if(line(1:23).eq.'_atom_site.label_seq_id') nresid=jtype
!       if(line(1:18).eq.'_atom_site.Cartn_x') nxcrd=jtype
!       if(line(1:18).eq.'_atom_site.Cartn_y') nycrd=jtype
!       if(line(1:18).eq.'_atom_site.Cartn_z') nzcrd=jtype
!       if(line(1:20).eq.'_atom_site.occupancy') noccp=jtype
!       if(line(1:23).eq.'_atom_site.label_alt_id') nalter=jtype
!       if(line(1:29).eq.'_atom_site.pdbx_PDB_model_num')nmrmodel = jtype

       if((line(1:20).eq.'_atom_site.group_PDB').and.
     1 (line(21:21).ne.'_')) natmhtm=jtype
       if((line(1:24).eq.'_atom_site.label_atom_id') .and.
     1 (line(25:25).ne.'_')) natm=jtype
       if((line(1:24).eq.'_atom_site.label_comp_id') .and.
     1 (line(25:25).ne.'_')) inres=jtype
       if((line(1:23).eq.'_atom_site.auth_asym_id') .and.
     1 (line(24:24).ne.'_')) nchn=jtype
       if((line(1:22).eq.'_atom_site.auth_seq_id') .and.
     1 (line(23:23).ne.'_')) nresid=jtype
       if((line(1:28).eq.'_atom_site.pdbx_PDB_ins_code') .and.
     1 (line(29:29).ne.'_')) nresid2=jtype
       if((line(1:18).eq.'_atom_site.Cartn_x') .and.
     1 (line(19:19).ne.'_')) nxcrd=jtype
       if((line(1:18).eq.'_atom_site.Cartn_y') .and.
     1 (line(19:19).ne.'_')) nycrd=jtype
       if((line(1:18).eq.'_atom_site.Cartn_z') .and.
     1 (line(19:19).ne.'_')) nzcrd=jtype
       if((line(1:20).eq.'_atom_site.occupancy') .and.
     1 (line(21:21).ne.'_')) noccp=jtype
       if((line(1:23).eq.'_atom_site.label_alt_id') .and.
     1 (line(24:24).ne.'_')) nalter=jtype
       if((line(1:29).eq.'_atom_site.pdbx_PDB_model_num') .and.
     1 (line(30:30).ne.'_')) nmrmodel=jtype
                read(4,1) line
                jtype=jtype+1
!--------------------------
          if(line(1:5).ne.'_atom') ik=jtype
                if(line(1:1).eq.'#') then
                  k=2
                  i=1
                  exit
                endif
                numatm=numatm+1
                if(numatm==ik-1)then
C               write(6,*) numatm,line
         if(natm.ne.0.and.inres.ne.0.and.nresid.ne.0.and.nxcrd.ne.0.
     1  .and.nycrd.ne.0.and.nzcrd.ne.0) then

          call stream(line,natm,atmnm,inres,res,nresid,idr,
     1   nchn,chain,nxcrd,xcoor,nycrd,ycoor,nzcrd,zcoor,noccp,occp,
     2   nmrmodel,natmhtm,static,nalter,alternate,modelfound,nresid2,
     3   inscode)
          ichain=index(chain,' ')
c          write(6,*) 'In READCIF',nresid,idr,nresid2,inscode,modelfound
          if(static.ne.1.and.(alternate.eq.'A'.or.alternate.eq.'.'))then
c          write(7,3) numatm,atmnm,res,chain,idr,xcoor,ycoor,zcoor,occp
            if(imolid.eq.1.or.molid(1:imolid).eq.chain(1:ichain)) then
              call hypothesis(atmnm(1:4),res(1:3),chain(1:4),idr,
     1 posei,xcoor,ycoor,zcoor,occp,npoorat,inscode)
            endif
c          write(7,3) numatm,atmnm,res,chain,idr,xcoor,ycoor,zcoor,occp
c          write(6,*) 'mresd',mresd
          endif
          endif
              endif
              enddo
         endif
         endif
         enddo
99      continue              
199     continue
1       format(a132)
3       format('ATOM',i7,2x,a3,1x,a3,1x,a4,i4,4x,3f8.3,f8.2)
!                    write(*,*) (mresd(ip),ip=70,80)
!        if(modelfound.eq.0) stop 9
        return
        end
!=====================================================================
!=====================================================================
        subroutine stream(data,nname,name,nresn,res,nresid,idr,nchn,
     1 chain,nx,xcr,ny,ycr,nz,zcr,no,occ,nmrmodel,natmhtm,static,nalter,
     2 alternate,modelfound,nresid2,inscode)
        integer c, co, modelfound
        character*132 data
        character*15 cda(30),name,res,chain,seq,atmhtm,alternate,nscode
        character*1 inscode
        integer hetatm,static
        common /options/cutang,cuteng,cutoff,hetatm,modelno
          c=2		! bhatta March 17, 2021
          co=1
          nd=1
          static=0
          call removedta(data)
          do while(c.le.132)
            if(data(c:c).eq.' '.and.data(c-1:c-1).ne.' ') then
              cda(nd)=data(co:c)
              c=c+1
              co=c
              nd=nd+1
            else
              c=c+1
            endif
          enddo
          if(natmhtm.ne.0) then
            atmhtm=cda(natmhtm)
            call removebl(atmhtm)
            alternate=cda(nalter)
            call removebl(alternate)
c        write(6,*) 'nalter',nalter,alternate
            if(atmhtm(1:6).eq.'HETATM'.and.hetatm.ne.1) then
              static=1
              return
            endif
          endif
          read(cda(nmrmodel),*) nmrmdl
          if(nmrmdl.eq.modelno) then   !!!!!   .and.(alternate.eq.'A'.or.alternate.eq.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     1   '.')) then
            modelfound=1
            name=cda(nname)
            call removebl(name)
            nrm=index(name,'"')
            if(nrm.gt.0) name(nrm:nrm)=' '  !  bhatta March 17, 2021
            res=cda(nresn)
            call removebl(res)
            nscode=cda(nresid2)
            call removebl(nscode)
            inscode=nscode
c            call removebl(inscode)
            read(cda(nresid),5,ERR=99)idr
c            read(cda(nresid),*)idr
c        write(6,*) data(1:55),' In Character:',cda(nresid),' in num',ird
            if(nchn.ne.0) chain=cda(nchn)
            call removebl(chain)
              read(cda(nx),*) xcr
              read(cda(ny),*) ycr
              read(cda(nz),*) zcr
          if(no.ne.0) read(cda(no),*) occ

        endif
1       format(a132)
2       format(a15)
3       format('ATOM',8x,a6,1x,a3,a3,i4,4x,3f8.3,f8.2)
4	format(a15)
5       format(I7)
99      continue
        return
        end
        subroutine removebl(variable)
        character*15 variable
        do while(variable(1:1).eq.' '.or.variable(1:1).eq.'"')
          variable(1:1)=variable(2:2)
          variable(2:2)=variable(3:3)
          variable(3:3)=variable(4:4)
          variable(4:4)=variable(5:5)
          variable(5:5)=variable(6:6)
          variable(6:6)=' '
        enddo
        return
        end

        subroutine removedta(variable)
        character*132 variable
        do while(variable(1:1).eq.' '.or.variable(1:1).eq.'"')
          variable(1:1)=variable(2:2)
          variable(2:2)=variable(3:3)
          variable(3:3)=variable(4:4)
          variable(4:4)=variable(5:5)
          variable(5:5)=variable(6:6)
          variable(6:6)=' '
        enddo
        return
        end

!==============================================
       subroutine hypothesis(atmnam,resn,chan,nores,posins,xa,
     1 ya,za,ocp,npoorat,inscode)

       character*1 posins,pos(900000),tmppos,inscode,tmpins,allins
       character*4 CHAIND,chan,pcd,chainid,pchaind,molid
       character*3 adevar,guavar,cytvar,uravar,resn,resid,RESD
	common/local/chainid
       character*4 atmnam,NAMED
       integer nores,atgc,presd,prd
       real xa,ya,za,ocp
       COMMON /NAMES/NAMED(900000),RESD(900000)
       COMMON /CHAINS/CHAIND(900000)
       common /occur/occ(900000)
       COMMON /ALPHA/XD(900000),YD(900000),ZD(900000)
       COMMON /BETA/IRESD(900000),ISTS(1000000),IENS(1000000),nres,
     1     allins(1000000),nforce,navar,ngvar,ncvar,nuvar,mresd(1000000)
       COMMON /GAMA/ATGC(1000000)
       COMMON /DELTA/l,lresno
       COMMON /DELTAD/tmppos,tmpins
       common /residue/resid(900000)
       common /num/kresd(1000000),prd(1000000)
       common /cum/pchaind(1000000),pcd(1000000),molid
       common /basenms/adevar(200),guavar(200),cytvar(200),uravar(200)
  
!------------------------------------
C       chainid='XXXX'
            do j=1,2
              if(resn(3:3).eq.' ') then
               resn(2:3)=resn(1:2)
               resn(1:1)=' '
              endif
            enddo
            nfound=0
C
C This may need to be modified to analyze model even if that is poorly
C defined
C
            if(ocp.eq.0.0e0.or.(ocp.ge.0.5.and.atmnam(4:4).ne.
     1                               'B'))then
               do i=1,navar
                  if(resn.eq.adevar(i)) then
                     resid(l) = '  A'
                     ATGC(l) = 1
                     nfound=1
                  endif
               enddo
               do i=1,ngvar
                  if(resn.eq.guavar(i)) then
                     resid(l) = '  G'
                     ATGC(l) = 2
                     nfound=1
                  endif
               enddo
               do i=1,ncvar
                  if(resn.eq.cytvar(i)) then
                     resid(l) = '  C'
                     ATGC(l) = 3
                     nfound=1
                  endif
               enddo
               do i=1,nuvar
                  if(resn.eq.uravar(i)) then
                     resid(l) = '  U'
                     ATGC(l) = 4
                     nfound=1
                  endif
               enddo
               if(resn.eq.'N6G') then
                     resid(l) = 'N6G'
                     ATGC(l) = 5
                     nfound=1
               endif
               if((resn.eq.'PSU').or.(resn.eq.'FHU')) then
                     resid(l)='PSU'
                     ATGC(l) = 6
                     nfound=1
               endif
               if(resn.eq.'QUO') then
                     resid(l)='QUO'
                     ATGC(l) = 7
                     nfound=1
               endif
       if((atmnam.eq.' C1''').and.(resid(l).ne.'  A').and.(resid(l)
     1.ne.'  G').and.(resid(l).ne.'  C').and.(resid(l).ne.'  U')
     1.and.(resid(l).ne.'PSU').and.(resid(l).ne.'N6G').and.(resid(l)
     1.ne.'QUO'))then
c               write(*,337)nores,resn,chan
               write(52,337)nores,resn,chan
               endif

337            format('#MODRES',3X,I5,2X,A3,2X,A4)
338    format('#='/'#HEADER==== Number of POORLY defined and REJECTED'
     1  'atoms:',i5' ====')
               if(nfound.eq.1)then
                 named(l)=atmnam
                 resd(l)=resn
                 chaind(l)=chan
                 iresd(l)=nores
                 pos(l)=posins
                 allins(l)=inscode
                 xd(l)=xa
                 yd(l)=ya
                 zd(l)=za
                 occ(l)=ocp
!------------------------------
c       write(6,*) atmnam,': ',resn,': chan:',chan,': ',nores,': ',
c     1 ':',l,': ',iresd(l),': ',lresno,': chainid::',
c     2 chainid,': ',chaind(l) !,':'
!	write(6,*) 'chainid:',chainid,'==chaind(l):',chaind(l)
!------------------------------
                 if(lresno.ne.iresd(l).or.tmppos.ne.pos(l).or.
     1   tmpins.ne.allins(l).or.chainid(1:4).ne.chaind(l)(1:4)) then  
                   nres=nres+1
                   ists(nres)=l
                   lresno=iresd(l)
                   tmppos=pos(l)
                   tmpins=allins(l)
                   chainid=chaind(l)
!	write(6,*) 'Updated chainid',chainid
                   kresd(nres)=lresno
                   mresd(nres)=lresno
                   pchaind(nres)=chainid
c              write(6,*) 'In Hypothesis',nres,lresno,kresd(nres),
c     1   ' &mRESD',mresd(nres)
!	write(6,*) 'chainid:',chainid,' lresno=',lresno,' nres=',nres
                 endif
                 iens(nres)=l
                 l = l+1
               endif
            else
              npoorat=npoorat+1
            endif
c       write(6,*) atmnam,':',resn,':',chan,':',nores,':',xa,
c     1 ya,za,ocp,':',l,':',iresd(l),':',lresno,':'
           return
           end
        subroutine renbase(bs,base)
        character*3 bs,adevar,guavar,cytvar,uravar
        character*1 base,allins
       COMMON /BETA/IRESD(900000),ISTS(1000000),IENS(1000000),nres,
     1   allins(1000000),nforce,navar,ngvar,ncvar,nuvar,mresd(1000000)
       common /basenms/adevar(200),guavar(200),cytvar(200),uravar(200)

   	do i=1,navar
	   if(bs.eq.adevar(i)) base='a'
	enddo
	do i=1,ngvar
	   if(bs.eq.guavar(i)) base='g'
	enddo
	do i=1,ncvar
	   if(bs.eq.cytvar(i)) base='c'
	enddo
	do i=1,nuvar
	   if(bs.eq.uravar(i)) base='u'
	enddo

        if(bs.eq.'  A') base='A'
        if(bs.eq.'  G') base='G'
        if(bs.eq.'  C') base='C'
        if(bs.eq.'  U') base='U'
        if(bs.eq.'PSU') base='u'
        if(bs.eq.'  T') base='T'
        return
        end
        function calcc1dis(i,j)
        character*4 named
        character*3 resd
        character*1 allins
       COMMON /NAMES/NAMED(900000),RESD(900000)
       COMMON /ALPHA/XD(900000),YD(900000),ZD(900000)
       COMMON /BETA/IRESD(900000),ISTS(1000000),IENS(1000000),nres,
     1   allins(1000000),nforce,navar,ngvar,ncvar,nuvar,mresd(1000000)
        do k=ists(i),iens(i)
          if(named(k).eq.'C1* ') then
            x11=xd(k)
            y11=yd(k)
            z11=zd(k)
          endif
        enddo
        do k=ists(j),iens(j)
          if(named(k).eq.'C1* ') then
            x12=xd(k)
            y12=yd(k)
            z12=zd(k)
          endif
        enddo
        calcc1dis=sqrt((x11-x12)**2 + (y11-y12)**2 + (z11-z12)**2)
        return
        end
