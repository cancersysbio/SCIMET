#!/bin/bash


WORKDIR='DIR'
TREEDIR='DIR/deme_trees'

cd $WORKDIR

NAME=`head -n $1 abc_parameters.txt | tail -n 1`
MUTRATE=`echo $NAME | awk '{print $1}'`
METTIME=`echo $NAME | awk '{print $2}'`
REPLICATE=`echo $NAME | awk '{print $3}'`


#echo $MUTRATE $METTIME $REPLICATE

source ~/.bashrc


./TumorMetTiming_abc_NN.py $MUTRATE $METTIME $TREEDIR/tumor3D_primary4samples_Nd${METTIME}.tree${REPLICATE} $TREEDIR/tumor3D_met4samples.tree${REPLICATE} $REPLICATE
./TumorMetTiming_abc_NS.py $MUTRATE $METTIME $TREEDIR/tumor3D_primary4samples_Nd${METTIME}.tree${REPLICATE} $TREEDIR/tumor3D_met4samples.tree${REPLICATE} $REPLICATE
./TumorMetTiming_abc_SN.py $MUTRATE $METTIME $TREEDIR/tumor3D_primary4samples_Nd${METTIME}.tree${REPLICATE} $TREEDIR/tumor3D_met4samples.tree${REPLICATE} $REPLICATE
./TumorMetTiming_abc_SS.py $MUTRATE $METTIME $TREEDIR/tumor3D_primary4samples_Nd${METTIME}.tree${REPLICATE} $TREEDIR/tumor3D_met4samples.tree${REPLICATE} $REPLICATE
