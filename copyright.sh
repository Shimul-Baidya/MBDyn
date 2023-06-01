#!/bin/sh
# $Header$

OLD=2021
NEW=`expr $OLD + 2`

for i in `grep -ErIl 'Copyright.*-'${OLD} .`; do 
	echo $i
	sed 's;\(Copyright.*\)'${OLD}';\1'${NEW}';' $i > x 
	mv -f x $i
done
