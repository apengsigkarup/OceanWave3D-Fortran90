#!/bin/sh
PRODUCT="OceanWave3D"
USER=apek
NUM=0.99.0
DIR="$PRODUCT"
svn export svn://$USER@svn.gbar.dtu.dk/$USER/OceanWave3D/tags/Release-$NUM ../$DIR-$NUM
tar czf $PRODUCT-$NUM.tar.gz ../$DIR-$NUM
mv $PRODUCT-$NUM.tar.gz ../.


