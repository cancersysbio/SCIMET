#!/bin/bash

WORKDIR='DIR'

cd $WORKDIR

NAME=`head -n $1 simulation_parameters.txt | tail -n 1`
MODEL=`echo $NAME | awk '{print $1}'`
REPLICATE=`echo $NAME | awk '{print $2}'`

#echo $MODEL $REPLICATE

source ~/.bashrc

./TumorSimul3D_LpLm_4models.py $REPLICATE $MODEL
