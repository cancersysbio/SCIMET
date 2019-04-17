#!/bin/bash


WORKDIR='DIR/abc'
TREEDIR='DIR/deme_trees'

cd $WORKDIR

NAME=`head -n $1 abc_parameters.txt | tail -n 1`
MUTRATE=`echo $NAME | awk '{print $1}'`
METTIME=`echo $NAME | awk '{print $2}'`
REPLICATE=`echo $NAME | awk '{print $3}'`


#echo $MUTRATE $METTIME $REPLICATE

source ~/.bashrc


./TumorSim3D_abc_NN.py $MUTRATE $METTIME $TREEDIR/tumor3D_primary${METTIME}.tree${REPLICATE} $TREEDIR/tumor3D_met4cores.tree${REPLICATE} $REPLICATE
./TumorSim3D_abc_NS.py $MUTRATE $METTIME $TREEDIR/tumor3D_primary${METTIME}.tree${REPLICATE} $TREEDIR/tumor3D_met4cores.tree${REPLICATE} $REPLICATE
./TumorSim3D_abc_SN.py $MUTRATE $METTIME $TREEDIR/tumor3D_primary${METTIME}.tree${REPLICATE} $TREEDIR/tumor3D_met4cores.tree${REPLICATE} $REPLICATE
./TumorSim3D_abc_SS.py $MUTRATE $METTIME $TREEDIR/tumor3D_primary${METTIME}.tree${REPLICATE} $TREEDIR/tumor3D_met4cores.tree${REPLICATE} $REPLICATE
