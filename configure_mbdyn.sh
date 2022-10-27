#!/bin/zsh
LD_LIBRARY_PATH="LD_LIBRARY_PATH=/usr/local/lib:/usr/lib64/"
LD_RUN_PATH="LD_RUN_PATH=/usr/local/lib"
PKG_CONFIG_PATH="PKG_CONFIG_PATH=/usr/local/lib/pkgconfig/"
LDFLAGS="LDFLAGS=\"-rdynamic -lsuitesparseconfig\""
CPPFLAGS="CPPFLAGS=\"-I/usr/include/suitesparse/\""

#CXXFLAGS="CXXFLAGS=\"-O0 -g\""
#FCFLAGS="FCFLAGS=\"-O0 -g\""
#CFLAGS="CFLAGS=\"-O0 -g\""


#CXXFLAGS="CXXFLAGS=\"-O2\""
#FCFLAGS="FCFLAGS=\"-O2\""
#CFLAGS="CFLAGS=\"-O2 \""

# added -fallow-argument-mismatch to cope with gcc10 and over
CXXFLAGS="CXXFLAGS=\"-O3 -w -g -std=c++11\""
FCFLAGS="FCFLAGS=\" -O3 -w -g -fallow-argument-mismatch\""
FFLAGS="FFLAGS=\"-O3 -w -g -fallow-argument-mismatch\""
CFLAGS="CFLAGS=\"-O3 -w -g -std=gnu17\""


MODULES="rotor_disc" #hid flightgear 
OPTIONS="--with-lapack --with-metis=no --enable-runtime-loading --with-blas=blas"

if [ -f "./configure" ]
then
	# aclocal
	# sed -e 's/\$wl--export-dynamic/\${wl}-export-dynamic/' aclocal.m4 > xxx && mv -f xxx aclocal.m4
# 	sed -e 's/\${wl}--export-dynamic/\${wl}-export-dynamic/' aclocal.m4 > xxx && mv -f xxx aclocal.m4
# 	libtoolize --automake --force --copy
# 	autoheader
# 	automake --foreign --add-missing --copy --force-missing
# 	autoconf 
	eval "$CFLAGS $FCFLAGS $FFLAGS $CXX $CXXFLAGS $CPPFLAGS $LD_LIBRARY_PATH \
		$LD_RUN_PATH $PKG_CONFIG_PATH $CPPFLAGS $CPPFLAGS \
		$LDFLAGS ./configure --with-module=\"$MODULES\" $OPTIONS"
else
	echo "CONFIGURE-MBDYN: could not locate configure script (check dir?). Aborting..."
	echo -e "\n\n"
	exit -1
fi
