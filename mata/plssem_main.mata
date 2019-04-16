*!plssem_main version 0.3.0
*!Written 15Jan2019
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

clear mata

set matastrict on
set matafavor speed

do plssem_classes.mata
do plssem_base.mata
do plssem_mga.mata
do plssem_rebus.mata
do plssem_fimix.mata
do plssem_gas.mata
do plssem_utils.mata

lmbuild lplssem.mlib, replace
