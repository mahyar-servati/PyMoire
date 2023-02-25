#!/bin/bash
#Author: Mahyar Servati
#Email: mahyar.servati@znu.ac.ir
#**********************************************************************************************************************
#Please add in bashrc or specify here Quantum_Espresso/bin, wannier and pseudo-potential paths
#bin=~/Quantum_ESPRESSO_6.6/qe-6.6/bin
#wandir=~/wannier90-3.1.0
ps=~/{pseudo potential directory}
#---------------------------------------------------------------------
#Choose prefix and folder name 
name=graphite
#---------------------------------------------------------------------
#Choose which calculations do you want. specify with y or n
scf=y
nscf=y
wannier=y
postwannier90=n
band=n
dos=n
#---------------------------------------------------------------------
#Do you want to remove save folder or not. specify with y or n
removesave=y
#---------------------------------------------------------------------
#If you want or not to run in parallel mode please specify with y or n, and write the number of cores you want to work.
parallel=y
core=16
#**********************************************************************************************************************
#**********************************************************************************************************************
#check failure function
check_failure () {
    if [ -z "$jobdone" ]
    then
        echo "encountered an error. please check output file."
    else
        echo "performed successfully."
    fi
}
#*********************************************************************
#Locate run files
PW="pw.x < "
BANDS="bands.x < "
PLOTBAND="plotband.x < "
DOS="dos.x < "
PROJWFC="projwfc.x < "
PW2WAN="pw2wannier90.x < "
WANPP="wannier90.x -pp "
WAN="wannier90.x "
POSTW90="postw90.x "
if test "$parallel" = "y" ; then
    PW="mpirun -np $core pw.x -inp "
    BANDS="mpirun -np $core bands.x -inp "
    DOS="mpirun -np $core dos.x -inp "
    PROJWFC="mpirun -np $core projwfc.x -inp "
    PW2WAN="mpirun -np $core pw2wannier90.x -inp "
    WAN="mpirun -np $core wannier90.x "
fi

if test "$bin" ; then
PW="$bin/pw.x < "
BANDS="$bin/bands.x < "
PLOTBAND="$bin/plotband.x < "
DOS="$bin/dos.x < "
PROJWFC="$bin/projwfc.x < "
PW2WAN="$bin/pw2wannier90.x < "
if test "$parallel" = "y" ; then
    PW="mpirun -np $core $bin/pw.x -inp "
    BANDS="mpirun -np $core $bin/bands.x -inp "
    DOS="mpirun -np $core $bin/dos.x -inp "
    PROJWFC="mpirun -np $core $bin/projwfc.x -inp "
    PW2WAN="mpirun -np $core $bin/pw2wannier90.x -inp "
fi
fi
if test "$wandir" ; then
WANPP="$wandir/wannier90.x -pp "
WAN="$wandir/wannier90.x "
POSTW90="$wandir/postw90.x "
if test "$parallel" = "y" ; then
    WAN="mpirun -np $core $wandir/wannier90.x "
fi
fi

#*********************************************************************
# self-consistent calculation
Firsttime=$(date +"%T")
cal_num=1
a=2.4599989204
mkdir out  #comment it in secound run
# displacing in a 60*60 mesh in a square space about 2*2 angstrom
for x in $(seq -30 1 30); do   #set displacing range in x axes  
dxFirsttime=$(date +"%T")
dx=`echo "scale=10; (1.1*$a/30.0)*$x" | bc -l` 

for y in $(seq -30 1 30); do   #set displacing range in y axes
dyFirsttime=$(date +"%T")
dy=`echo "scale=10; (1.1*$a/30.0)*$y" | bc -l`
echo "-------------------------calculation number $cal_num begins-------------------------"

foldername="$name-$x-$y"        #set folder name pattern- my first folder is (graphite-30--30)
if [ ! -d $foldername ] ; then
    mkdir $foldername
fi
cd $foldername

if test "$scf" = "y" ; then
    cat > $name.scf << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = '$ps',
    outdir='./'
    prefix='$name'
 /
 &system
                   vdw_corr = 'grimme-d2'
                       ibrav = 0,
                   celldm(1) = 1.8897261249935,
                         nat = 4
                        ntyp = 1,

    occupations='smearing', smearing='cold', degauss=0.02
    ecutwfc = 60
 /
 &electrons
    startingwfc='random'
    diagonalization='cg'
    conv_thr = 1.0e-8
 /
CELL_PARAMETERS cubic
    2.1304215583   -1.2299994602    0.0000000000
    0.0000000000    2.4599989204    0.0000000000
    0.0000000000    0.0000000000    12.000000000
ATOMIC_SPECIES
C  12.011    C.pbe-n-rrkjus_psl-1-0-0.UPF
ATOMIC_POSITIONS alat
     C    `echo "scale=2; 0.0000000000+$dx" | bc -l`    `echo "scale=2; 0.0000000000+$dy" | bc -l`    7.625000000
     C    `echo "scale=2; 0.7101405194+$dx" | bc -l`    `echo "scale=2; 1.2299994602+$dy" | bc -l`    7.625000000
     C    0.0000000000    0.0000000000    4.3750000000
     C    1.4202810389    0.0000000000    4.3750000000
K_POINTS (automatic)
10 10 1 0 0 0


EOF
    echo "running the scf calculation for $x$y ..."
    `$PW $name.scf > $name.scf.out`
    jobdone=`awk '/JOB DONE./ { print $0 }' $name.scf.out`
    echo -n "Step 1: scf calculation "
    check_failure

    TE=`awk '/!/ { print $5 }' $name.scf.out`
    echo "Total energy is $TE"
    echo "$dx $dy $TE" >> ../out/Total_Energy.dat

    HOMO=`awk '/highest occupied/ { print $7 }' $name.scf.out`
    LUMO=`awk '/highest occupied/ { print $8 }' $name.scf.out`
    echo "$dx $dy $HOMO $LUMO" >> ../out/HOMO_LUMO.dat

    Evdw=`awk '/Dispersion Correction/ { print $4 }' $name.scf.out`
    echo "$dx $dy $Evdw" >> ../out/VdW_Energy.dat

    xc=`awk '/xc contribution/ { print $4 }' $name.scf.out`
    echo "$dx $dy $xc" >> ../out/xc_Energy.dat
fi
#*******************************************************************
# nscf calculation
if test "$nscf" = "y" ; then
    cat > $name.nscf << EOF
 &control
    calculation='nscf'
    pseudo_dir = '$ps',
    outdir='./'
    prefix='$name'
 /
 &system
                   vdw_corr = 'grimme-d2'
                       ibrav = 0,
                   celldm(1) = 1.8897261249935,
                         nat = 4,
                        ntyp = 1,
    nbnd=20,nosym=.true.,
    occupations='smearing', smearing='cold', degauss=0.02
    ecutwfc = 60
 /
 &electrons
    startingwfc='random'
    diagonalization='cg'
    conv_thr = 1.0e-8
    mixing_beta = 0.2
 /
CELL_PARAMETERS cubic
    2.1304215583   -1.2299994602    0.0000000000
    0.0000000000    2.4599989204    0.0000000000
    0.0000000000    0.0000000000    12.000000000
ATOMIC_SPECIES
C  12.011    C.pbe-n-rrkjus_psl-1-0-0.UPF
ATOMIC_POSITIONS alat
     C    `echo "scale=2; 0.0000000000+$dx" | bc -l`    `echo "scale=2; 0.0000000000+$dy" | bc -l`    7.625000000
     C    `echo "scale=2; 0.7101405194+$dx" | bc -l`    `echo "scale=2; 1.2299994602+$dy" | bc -l`    7.625000000
     C    0.0000000000    0.0000000000    4.3750000000
     C    1.4202810389    0.0000000000    4.3750000000
K_POINTS crystal
36
  0.00000000  0.00000000  0.00000000  2.777778e-02
  0.00000000  0.16666667  0.00000000  2.777778e-02
  0.00000000  0.33333333  0.00000000  2.777778e-02
  0.00000000  0.50000000  0.00000000  2.777778e-02
  0.00000000  0.66666667  0.00000000  2.777778e-02
  0.00000000  0.83333333  0.00000000  2.777778e-02
  0.16666667  0.00000000  0.00000000  2.777778e-02
  0.16666667  0.16666667  0.00000000  2.777778e-02
  0.16666667  0.33333333  0.00000000  2.777778e-02
  0.16666667  0.50000000  0.00000000  2.777778e-02
  0.16666667  0.66666667  0.00000000  2.777778e-02
  0.16666667  0.83333333  0.00000000  2.777778e-02
  0.33333333  0.00000000  0.00000000  2.777778e-02
  0.33333333  0.16666667  0.00000000  2.777778e-02
  0.33333333  0.33333333  0.00000000  2.777778e-02
  0.33333333  0.50000000  0.00000000  2.777778e-02
  0.33333333  0.66666667  0.00000000  2.777778e-02
  0.33333333  0.83333333  0.00000000  2.777778e-02
  0.50000000  0.00000000  0.00000000  2.777778e-02
  0.50000000  0.16666667  0.00000000  2.777778e-02
  0.50000000  0.33333333  0.00000000  2.777778e-02
  0.50000000  0.50000000  0.00000000  2.777778e-02
  0.50000000  0.66666667  0.00000000  2.777778e-02
  0.50000000  0.83333333  0.00000000  2.777778e-02
  0.66666667  0.00000000  0.00000000  2.777778e-02
  0.66666667  0.16666667  0.00000000  2.777778e-02
  0.66666667  0.33333333  0.00000000  2.777778e-02
  0.66666667  0.50000000  0.00000000  2.777778e-02
  0.66666667  0.66666667  0.00000000  2.777778e-02
  0.66666667  0.83333333  0.00000000  2.777778e-02
  0.83333333  0.00000000  0.00000000  2.777778e-02
  0.83333333  0.16666667  0.00000000  2.777778e-02
  0.83333333  0.33333333  0.00000000  2.777778e-02
  0.83333333  0.50000000  0.00000000  2.777778e-02
  0.83333333  0.66666667  0.00000000  2.777778e-02
  0.83333333  0.83333333  0.00000000  2.777778e-02  
EOF
    echo "running the nscf calculation for $x$y ..."
    `$PW $name.nscf > $name.nscf.out`
    jobdone=`awk '/JOB DONE./ { print $0 }' $name.nscf.out`
    echo -n "Step 2: nscf calculation "
    check_failure
Ef=`awk '/the Fermi energy is/ { print $5 }' $name.nscf.out`
echo "Fermi energy is $Ef"
echo "$dx $dy $Ef" >> ../out/Fermi_Energy.dat
fi

Ef=`awk '/Fermi/ { print $5 }' $name.nscf.out`
#********************************************************************
#nnkp calculation
if test "$wannier" = "y" ; then
    cat > $name.win <<EOF
num_bands         =   20       
num_wann          =   10 

dis_win_max       = 19.2d0
dis_froz_max      =  5.0d0
dis_num_iter      =  300
num_iter          = 200
num_print_cycles  = 10

guiding_centres = true
write_hr = true
write_xyz = True
translate_home_cell = False
bands_plot = true

!spinors = true

begin kpoint_path
G 0.0000000000     0.0000000000     0.0000000000  M 0.5000000000     0.0000000000     0.0000000000
M 0.5000000000     0.0000000000     0.0000000000  K 0.3333333333     0.3333333333     0.0000000000 
K 0.3333333333     0.3333333333     0.0000000000  G 0.0000000000     0.0000000000     0.0000000000
G 0.0000000000     0.0000000000     0.0000000000  A 0.0000000000     0.0000000000     0.5000000000

end kpoint_path


Begin Atoms_Cart
     C1    `echo "scale=2; 0.0000000000+$dx" | bc -l`    `echo "scale=2; 0.0000000000+$dy" | bc -l`    7.625000000
     C3    `echo "scale=2; 0.7101405194+$dx" | bc -l`    `echo "scale=2; 1.2299994602+$dy" | bc -l`    7.625000000
     C2    0.0000000000    0.0000000000    4.3750000000
     C4    1.4202810389    0.0000000000    4.3750000000
End Atoms_Cart
    
Begin Projections     
C1:sp2;pz
C3:pz
C2:sp2;pz
C4:pz
End Projections       
    
Begin Unit_Cell_Cart
    2.1304215583   -1.2299994602    0.0000000000
    0.0000000000    2.4599989204    0.0000000000
    0.0000000000    0.0000000000    12.000000000
End Unit_Cell_Cart


mp_grid      = 6 6 1  


begin kpoints
  0.00000000  0.00000000  0.00000000
  0.00000000  0.16666667  0.00000000
  0.00000000  0.33333333  0.00000000
  0.00000000  0.50000000  0.00000000
  0.00000000  0.66666667  0.00000000
  0.00000000  0.83333333  0.00000000
  0.16666667  0.00000000  0.00000000
  0.16666667  0.16666667  0.00000000
  0.16666667  0.33333333  0.00000000
  0.16666667  0.50000000  0.00000000
  0.16666667  0.66666667  0.00000000
  0.16666667  0.83333333  0.00000000
  0.33333333  0.00000000  0.00000000
  0.33333333  0.16666667  0.00000000
  0.33333333  0.33333333  0.00000000
  0.33333333  0.50000000  0.00000000
  0.33333333  0.66666667  0.00000000
  0.33333333  0.83333333  0.00000000
  0.50000000  0.00000000  0.00000000
  0.50000000  0.16666667  0.00000000
  0.50000000  0.33333333  0.00000000
  0.50000000  0.50000000  0.00000000
  0.50000000  0.66666667  0.00000000
  0.50000000  0.83333333  0.00000000
  0.66666667  0.00000000  0.00000000
  0.66666667  0.16666667  0.00000000
  0.66666667  0.33333333  0.00000000
  0.66666667  0.50000000  0.00000000
  0.66666667  0.66666667  0.00000000
  0.66666667  0.83333333  0.00000000
  0.83333333  0.00000000  0.00000000
  0.83333333  0.16666667  0.00000000
  0.83333333  0.33333333  0.00000000
  0.83333333  0.50000000  0.00000000
  0.83333333  0.66666667  0.00000000
  0.83333333  0.83333333  0.00000000
End Kpoints
EOF

    echo "running the wannier estimate for $x$y ..."
    `$WANPP $name`
    nnkp=`awk '/nnkp written/ { print $0 }' $name.wout`
    if [ -z "$nnkp" ]
    then
        echo "Step 3: nnkp calculation encountered an error. please check output file."
    else
        echo "Step 3: nnkp calculation performed successfully."
    fi
#*******************************************************************
#pw2wan calculation
    cat > $name.pw2wan <<EOF
&inputpp
  prefix='$name'
  spin_component = 'none'
  write_mmn = .true.
  write_amn = .true.
!  write_spn = .true.
!  write_unk = .false.
  wan_mode='standalone'
  seedname = '$name'
/
EOF
    echo "running the pw2wan transformation calculation for $x$y ..."
    `$PW2WAN $name.pw2wan > $name.pw2wan.out`
    jobdone=`awk '/JOB DONE./ { print $0 }' $name.pw2wan.out`
    echo -n "Step 4: pw2wan calculation "
    check_failure

#********************************************************************
# wannier90 calculation
    echo "running the wannier optimize calculation for $x$y ..."
    `$WAN $name`
    alldone=`awk '/All done:/ { print $0 }' $name.wout`
    if [ -z "$alldone" ]
    then
        echo "Step 5: wannier optimize calculation encountered an error. please check output file."
    else
        echo "Step 5: wannier optimize calculation performed successfully."
    fi
	bandnum=`awk 'NR == 2 {print $1}' "$name"_hr.dat`	
	wsnum=`awk 'NR == 3 {print $1}' "$name"_hr.dat`
	skipline=`echo "scale=0; ($wsnum/15)+4" | bc -l`
	for l1 in $(seq 1 1 $bandnum); do
	for l2 in $(seq 1 1 $bandnum); do
	hrmat="$l1-$l2.dat"
	hop=`awk -v var="$skipline" -v l1="$l1" -v l2="$l2" 'NR > var {if($4 == l1 && $5 == l2) print $6 }' "$name"_hr.dat`
	echo "$dx $dy" $(echo "${hop[@]}") >> ../out/$hrmat
	done
	done
    echo "tight-binding parameters extracted from hr file in term of coupling orbitals and saved in out folder"
fi
#********************************************************************
# postw90 calculation
if test "$postwannier90" = "y" ; then
    echo "running the postw90 calculation for $foldername ..."
    `$POSTW90 $name`
    alldone=`awk '/All done:/ { print $0 }' $name.wpout`
    if [ -z "$alldone" ]
    then
        echo "Step 5: postw90 calculation encountered an error. please check output file."
    else
        echo "Step 5: postw90 calculation performed successfully."
    fi
fi
#*******************************************************************
# bands calculation
if test "$band" = "y" ; then
    cat > $name.bands << EOF
 &control
    calculation='bands'
    pseudo_dir = '$ps',
    outdir='./'
    prefix='$name'
 /
 &system
                   vdw_corr = 'grimme-d2'
                       ibrav = 0,
                   celldm(1) = 1.8897261249935,
                         nat = 4,
                        ntyp = 1,
    nbnd=20,nosym=.true.,
    occupations='smearing', smearing='cold', degauss=0.02
    ecutwfc = 60
 /
 &electrons
    startingwfc='random'
    diagonalization='cg'
    conv_thr = 1.0e-8
    mixing_beta = 0.2
 /
CELL_PARAMETERS cubic
    2.1304215583   -1.2299994602    0.0000000000
    0.0000000000    2.4599989204    0.0000000000
    0.0000000000    0.0000000000    12.000000000
ATOMIC_SPECIES
C  12.011    C.pbe-n-rrkjus_psl-1-0-0.UPF
ATOMIC_POSITIONS alat
     C    `echo "scale=2; 0.0000000000+$dx" | bc -l`    `echo "scale=2; 0.0000000000+$dy" | bc -l`    7.625000000
     C    `echo "scale=2; 0.7101405194+$dx" | bc -l`    `echo "scale=2; 1.2299994602+$dy" | bc -l`    7.625000000
     C    0.0000000000    0.0000000000    4.3750000000
     C    1.4202810389    0.0000000000    4.3750000000
K_POINTS {crystal_b}
5
0.000000000  0.000000000  0.000000000 101 !Gamma
0.500000000  0.000000000  0.000000000 58  !M
0.333333333  0.333333333  0.000000000 115 !K
0.000000000  0.000000000  0.000000000 31  !Gamma
0.000000000  0.000000000  0.500000000 0   !A
EOF
    echo "running the bands calculation for $x$y ..."
    `$PW $name.bands > $name.bands.out`

    jobdone=`awk '/JOB DONE./ { print $0 }' $name.bands.out`
    echo -n "Step 6: bands calculation "
    check_failure

#*******************************************************************
# ppbands calculation

    cat > $name.ppbands << EOF
&BANDS
filband = '`echo $name`_bandsx.dat'
!lsym = .TRUE.
outdir = './'
prefix = '$name'
/
EOF
    echo "running the ppbands calculation for $x$y ..."
    `$BANDS $name.ppbands > $name.ppbands.out`

    jobdone=`awk '/JOB DONE./ { print $0 }' $name.ppbands.out`
    echo -n "Step 7: ppbands calculation "
    check_failure
#*******************************************************************
# plotband calculation
#    cat > $name.plotband << EOF
#`echo $name`_bandsx.dat
#-25 15
#$name.bands.xmgr
#$name.bands.ps
#$Ef
#2 $Ef
#EOF
#    echo "running the plotband calculation for $dx$dy ..."
#    `$PLOTBAND $name.plotband > $name.plotband.out`

#    echo "Step 8: plotband calculation performed successfully.";
fi

#*******************************************************************
# dos calculation
if test "$dos" = "y" ; then
    cat > $name.dos << EOF
&dos
    outdir='./'
    prefix='$name'
    fildos='$name.dos.result',
    Emin=-20.0, Emax=15.0, DeltaE=0.1
 /
EOF
    echo "running the dos calculation for $x$y ..."
    `$DOS $name.dos > $name.dos.out`
    jobdone=`awk '/JOB DONE./ { print $0 }' $name.dos.out`
    echo -n "Step 9: dos calculation "
    check_failure

#*******************************************************************
# projwfc calculation

    cat > $name.projwfc << EOF
&projwfc
    outdir='./'
    prefix='$name'
    Emin=-20.0, Emax=15.0, DeltaE=0.1
    ngauss=1, degauss=0.02
 /

EOF
    echo "running the projwfc calculation for $x$y ..."
    `$PROJWFC $name.projwfc > $name.projwfc.out`
    jobdone=`awk '/JOB DONE./ { print $0 }' $name.projwfc.out`
    echo -n "Step 10: projwfc calculation "
    check_failure
fi
dyLasttime=$(date +"%T")
echo "************ $x$y calculations done between $dyFirsttime and $dyLasttime ************"
if test "$removesave" = "y" ; then
rm -r $name.save
fi
cd ../
cal_num=$(($cal_num +1))
done          #dy done
dxLasttime=$(date +"%T")
echo "************ All dx= $x calculations done between $dxFirsttime and $dxLasttime ************"
done          #dx done
Totaltime=$(date +"%T")
echo "************ All calculations done between $Firsttime and $Totaltime ************"
