#!/bin/csh

rm -f Makefile >& /dev/null

######
###### QUESTION #1
######
echo "################################################################################"
echo "Do you want a diatomic (diat) fit, a triatomic (triat) or tetratomic (tetra) ?"
echo "################################################################################"

set MOL=$<

if ( $MOL == "diat" ) then 

#  echo "diatomic"

else if ( $MOL == "triat" ) then

#  echo "triatomic"

else if ( $MOL == "tetra" ) then

#  echo "tetratomic"

else 
 
  echo "Please, provide a valid parameter: "diat" for diatomics, "triat" for triatomics or "tetra" for tetratomics"
  exit

endif
######
###### QUESTION #2
######
echo "################################################################################"
echo "gfortran is used by default. Do you agree ? [y/n]"  
echo "################################################################################"

set COMPTYP=$<

if ( $COMPTYP == "y" ) then 

  set FC="gfortran"
  set CFLAGS="-g -O3" 
#
# set CFLAGS="-g -O3 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace" # for debug
#
  set LDFLAGS="-g" 

else if ( $COMPTYP == "n" ) then

  set FC="XXX"
  set CFLAGS="XXX"
  set LDFLAGS="XXX" 

  echo "Please, configure latter the Makefile with your compiler options"

else
 
  echo "Please, provide a proper answer: yes (y) or no (n)"
  exit

endif
######
###### QUESTION #3
######
echo "################################################################################"
echo "Is a basis set calibration (bas) planned, a polynomial fit (pol) or a direct fit (dirfit)?"
echo "################################################################################"

set JOB=$<

if ( $JOB == "bas" ) then 
 
else if ( $JOB == "pol" ) then

else if ( $JOB == "dirfit" ) then

else
 
  echo "Please, provide a proper answer: "bas" for basis set calibration, "pol" for polynomial fit or "dirfit" for direct fit"
  exit

endif

######
###### QUESTION #4
######
echo "################################################################################"
echo "Are additional subroutines to be included into the program ? [y/n]"
echo "################################################################################"

set ADDSUB=$< 

if ( $ADDSUB == "y" ) then 

  echo "################################################################################"
  echo "If yes, specify its name, e.g., XXX.f or XXX.f90"
  echo "################################################################################"

  set ADDSUBNAME=$<

  if ($ADDSUBNAME:e == "f" || $ADDSUBNAME:e == "f90" || $ADDSUBNAME:e == "c") then 

  else
  echo "Please, provide a valid name for the additional module, e.g., XXX.f or XXX.f90"
  exit
  endif

  if (-e $ADDSUBNAME) then
  else
  echo "Please, check if this module is in this directory" 
  exit

  exit

  endif
 
else if ( $ADDSUB == "n" ) then

  set ADDSUBNAME=" "

else
 
  echo "Please, provide a proper answer: yes (y) or no (n)"
  exit

endif

######
###### LOCATING MODULES
######

set i = 0

set ARGVFILE = (`seq 1 1 23`)
set DIRFILE = (`seq 1 1 23`)

foreach FILE (MODULE_COMMON_VAR.f90 MAIN.f90 READING.f90 TOFIT.f90 LSFIT.f90 LMDIF.f CHIPR_DIAT.f90 CHIPR_TRIAT.f90 CHIPR_TETRA.f90 BASIS_SET.f90 WEIGHT_FUNC.f90 SUM2BD.f90 SUM23BD.f90 FINAL_ANALYSIS.f90 PROP_DIAT.f90 PROP_TRIAT.f90 PROP_TETRA.f90 STRAT_RMSD.f90 SA.f HARM_VIB_SA.f90 DIRECT_FIT.f90 LEVEL.f $ADDSUBNAME)

@ i = $i + 1

set ARGVFILE[$i] = $FILE

#echo "$ARGVFILE[$i]"

if ( -e $ARGVFILE[$i] ) then 

  echo "################################################################################"
  echo "$ARGVFILE[$i] will be set as the one contained in this directory"
  echo "################################################################################"

  set DIRFILE[$i] = "./$ARGVFILE[$i]:r"
#  echo "$DIRFILE[$i]"

else if ( ! -e $ARGVFILE[$i] ) then 

  set DIRFILE[$i] = "../CHIPR-4.0_SOURCE_CODE/$ARGVFILE[$i]:r"
#  echo "$DIRFILE[$i]"

endif

end

echo ""

if ($JOB == "bas" && $MOL == "diat") then

   echo "################################################################################"
   echo "Please, make sure that you have in this directory the following files"
   echo "################################################################################"
   echo "1). coeffs_guess.txt"
   echo "2). PLOT_BASIS_OPT_CHIPR_DIAT.gnu"
#   echo "3). Any additional subroutine, e.g., XXX.f90, you want to include or replace locally with any one contained in ../#CHIPR-4.0_SOURCE_CODE"

else if ($JOB == "pol" && $MOL == "diat") then

   echo "################################################################################"
   echo "Please, make sure that you have in this directory the following files"
   echo "################################################################################"
   echo "1). coeffs_guess.txt"
   echo "2). abinitio_data.txt"
   echo "3). mol_gm_guess.txt"
   echo "4). PLOT_POL_OPT_CHIPR_DIATOM.gnu"
   echo "5). "WEIGHT_FUNC.f90" if you want to set weights automatically"
#   echo "6). Any additional subroutine, e.g., XXX.f90, you want to include or replace locally with any one contained in ../#CHIPR-4.0_SOURCE_CODE"

else if ($JOB == "dirfit" && $MOL == "diat") then

   echo "################################################################################"
   echo "Please, make sure that you have in this directory the following files"
   echo "################################################################################"
   echo "1). coeffs_guess.txt"
   echo "2). abinitio_data.txt"
   echo "3). mol_gm_guess.txt"
   echo "4). spec_data.txt"
   echo "5). PLOT_POL_OPT_CHIPR_DIATOM.gnu"
   echo "6). "WEIGHT_FUNC.f90" if you want to set weights automatically"
#   echo "7). Any additional subroutine, e.g., XXX.f90, you want to include or replace locally with any one contained in ../#CHIPR-4.0_SOURCE_CODE"


else if ($JOB == "bas" && $MOL == "triat") then

   echo "################################################################################"
   echo "Please, make sure that you have in this directory the following files"
   echo "################################################################################"
   echo "1). coeffs_guess.txt"
   echo "2). PLOT_BASIS_OPT_CHIPR_TRIAT.gnu"
   echo "3). "SUM2BD.f90""
   echo "4). Any additional subroutine, e.g., XXX.f90, you want to include or replace locally with any one contained in ../CHIPR-4.0_SOURCE_CODE"

else if ($JOB == "pol" && $MOL == "triat") then

   echo "################################################################################"
   echo "Please, make sure that you have in this directory the following files"
   echo "################################################################################"
   echo "1). coeffs_guess.txt"
   echo "2). abinitio_data.txt"
   echo "3). mol_gm_guess.txt"
   echo "4). PLOT_POL_OPT_CHIPR_TRIATOM.gnu"
   echo "5). "SUM2BD.f90""
   echo "6). "WEIGHT_FUNC.f90" if you want to set weights automatically"
   echo "7). Any additional subroutine, e.g., XXX.f90, you want to include or replace locally with any one contained in ../CHIPR-4.0_SOURCE_CODE"

else if ($JOB == "bas" && $MOL == "tetra") then

   echo "################################################################################"
   echo "Please, make sure that you have in this directory the following files"
   echo "################################################################################"
   echo "1). coeffs_guess.txt"
   echo "2). PLOT_BASIS_OPT_CHIPR_TETRA.gnu"
   echo "3). "SUM23BD.f90""
   echo "4). Any additional subroutine, e.g., XXX.f90, you want to include or replace locally with any one contained in ../CHIPR-4.0_SOURCE_CODE"

else if ($JOB == "pol" && $MOL == "tetra") then

   echo "################################################################################"
   echo "Make sure that you have in this directory the following files"
   echo "################################################################################"
   echo "1). coeffs_guess.txt"
   echo "2). abinitio_data.txt"
   echo "3). mol_gm_guess.txt"
   echo "4). PLOT_POL_OPT_CHIPR_TETRA.gnu"
   echo "5). "SUM23BD.f90""
   echo "6). "WEIGHT_FUNC.f90" if you want to set weights automatically"
   echo "7). Any additional subroutine, e.g., XXX.f90, you want to include or replace locally with any one contained in ../CHIPR-4.0_SOURCE_CODE"

endif

echo ""

   echo "################################################################################"
   echo "Read to compile..."
   echo "################################################################################"
   echo "First, take a look in the Makefile so generated." 
   echo "If it is OK, then proceed to the compilation by typing the commands:"
   echo "1). make"
   echo "2). make clean"
   echo "################################################################################"
   echo "################################################################################"

cat << EOF > Makefile
FC=$FC
CFLAGS=$CFLAGS
LDFLAGS=$LDFLAGS

OBJS=   $ARGVFILE[1]:r.o \
	$ARGVFILE[2]:r.o \
	$ARGVFILE[3]:r.o \
	$ARGVFILE[4]:r.o \
	$ARGVFILE[5]:r.o \
	$ARGVFILE[6]:r.o \
	$ARGVFILE[7]:r.o \
	$ARGVFILE[8]:r.o \
	$ARGVFILE[9]:r.o \
	$ARGVFILE[10]:r.o \
	$ARGVFILE[11]:r.o \
	$ARGVFILE[12]:r.o \
	$ARGVFILE[13]:r.o \
	$ARGVFILE[14]:r.o \
	$ARGVFILE[15]:r.o \
	$ARGVFILE[16]:r.o \
	$ARGVFILE[17]:r.o \
	$ARGVFILE[18]:r.o \
	$ARGVFILE[19]:r.o \
	$ARGVFILE[20]:r.o \
	$ARGVFILE[21]:r.o \
	$ARGVFILE[22]:r.o \
EOF

if ( $ADDSUB == "y" ) then 

  echo "	$ARGVFILE[23]" >> Makefile
 
endif
	

cat << EOF >> Makefile

CHIPR.x : \$(OBJS) 
	\$(FC) \$(LDFLAGS) \$(OBJS) -o \$@

MODULE_COMMON_VAR.o : $DIRFILE[1].f90 
	\$(FC) \$(CFLAGS) -c $DIRFILE[1].f90

MAIN.o : $DIRFILE[2].f90 
	\$(FC) \$(CFLAGS) -c $DIRFILE[2].f90 

READING.o : $DIRFILE[3].f90 
	\$(FC) \$(CFLAGS) -c $DIRFILE[3].f90 

TOFIT.o :: $DIRFILE[4].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[4].f90

LSFIT.o :: $DIRFILE[5].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[5].f90

LMDIF.o : $DIRFILE[6].f
	\$(FC) \$(CFLAGS) -c $DIRFILE[6].f

CHIPR_DIAT.o :: $DIRFILE[7].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[7].f90

CHIPR_TRIAT.o :: $DIRFILE[8].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[8].f90

CHIPR_TETRA.o :: $DIRFILE[9].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[9].f90

BASIS_SET.o : $DIRFILE[10].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[10].f90

WEIGHT_FUNC.o : $DIRFILE[11].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[11].f90

SUM2BD.o :: $DIRFILE[12].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[12].f90

SUM23BD.o :: $DIRFILE[13].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[13].f90

FINAL_ANALYSIS.o :: $DIRFILE[14].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[14].f90

PROP_DIAT.o : $DIRFILE[15].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[15].f90

PROP_TRIAT.o : $DIRFILE[16].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[16].f90

PROP_TETRA.o : $DIRFILE[17].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[17].f90

STRAT_RMSD.o : $DIRFILE[18].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[18].f90

SA.o : $DIRFILE[19].f
	\$(FC) \$(CFLAGS) -c $DIRFILE[19].f

HARM_VIB_SA.o : $DIRFILE[20].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[20].f90

DIRECT_FIT.o : $DIRFILE[21].f90
	\$(FC) \$(CFLAGS) -c $DIRFILE[21].f90

LEVEL.o : $DIRFILE[22].f
	\$(FC) \$(CFLAGS) -c $DIRFILE[22].f

clean :
	rm $ARGVFILE[1]:r.o \
	$ARGVFILE[2]:r.o \
	$ARGVFILE[3]:r.o \
	$ARGVFILE[4]:r.o \
	$ARGVFILE[5]:r.o \
	$ARGVFILE[6]:r.o \
	$ARGVFILE[7]:r.o \
	$ARGVFILE[8]:r.o \
	$ARGVFILE[9]:r.o \
	$ARGVFILE[10]:r.o \
	$ARGVFILE[11]:r.o \
	$ARGVFILE[12]:r.o \
	$ARGVFILE[13]:r.o \
	$ARGVFILE[14]:r.o \
	$ARGVFILE[15]:r.o \
	$ARGVFILE[16]:r.o \
	$ARGVFILE[17]:r.o \
	$ARGVFILE[18]:r.o \
	$ARGVFILE[19]:r.o \
	$ARGVFILE[20]:r.o \
	$ARGVFILE[21]:r.o \
	$ARGVFILE[22]:r.o

EOF

exit 
