#!/usr/bin/env bash

mkdir -p $1/{cm,current/{calibration/inc_split,flow/{BL,VX},state/{BL,VX}},metadata,params/{seeds,sets},raw_output/traj}
touch $1/{cm,current/{calibration/inc_split,flow/{BL,VX},state/{BL,VX}},metadata,params/{seeds,sets},raw_output/traj}/.mkfolder
