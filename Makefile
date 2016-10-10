#
##  Makefile for mf_exciton, Carlo Motta, 2015
#
##  This Makefile works with GNU make. NOTE that on some proprietary platforms,
#  different versions of the ancient Unix make may not work - in that case, please
#  #  install GNU make - if it's not already there as "gmake".
#
##  MANDATORY REQUIREMENTS:
#   the only requirement is the Lapack library
#
SOURCECODE = mf_exciton_3d.f90
all:
	gfortran  ${SOURCECODE} -o exe /usr/lib/lapack/liblapack.so.3
