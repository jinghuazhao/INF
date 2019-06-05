module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.5.0-ICC-MKL
tar xvfz gap_1.2.1.tar.gz
cd gap/src
gcc -I/services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/include -L/services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib/release -L/services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker /services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib/release -Xlinker -rpath -Xlinker /services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib -Xlinker -rpath -Xlinker /opt/intel/mpi-rt/2017.0.0/intel64/lib/release -Xlinker -rpath -Xlinker /opt/intel/mpi-rt/2017.0.0/intel64/lib -lmpifort -lmpi -ldl -lrt -lpthread -L/services/tools/intel/perflibs/2019//compilers_and_libraries_2019.0.117/linux/mpi/intel64/libfabric/lib -fPIC -c *.c *.f
gcc -shared -L/services/tools/R/3.5.0-ICC-MKL/lib64/R/lib -L/usr/local/lib64 -o gap.so 2k.o 2ld.o cline.o gcontrol_c.o gcx.o gif_c.o hap_c.o hwe.hardy.o kin.morgan.o makeped_c.o mia.o muvar.o package_native_routine_registration_skeleton.o pfc.o pfc.sim.o pgc_c.o whscore.o -L/usr/lib/gcc/x86_64-redhat-linux/4.8.2 -lgfortran -lm -lquadmath -L/services/tools/R/3.5.0-ICC-MKL/lib64/R/lib -lR
cd -
R CMD INSTALL gap
