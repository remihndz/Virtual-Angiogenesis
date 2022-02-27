#!/bin/bash

# Runs the same simulation with three different 'minimum bifurcation angle' parameter
rm LogFiles/*
for P in {1..10} ;
do
    mkdir ../perfused_geometries/sim${P} 
    time ./dcco_svp ConfigFiles/Baseline.conf > LogFiles/sim${P}.log 
    mv ../perfused_geometries/Baseline.* ../perfused_geometries/sim${P}/ 
    cp LogFiles/sim${P}.log ../perfused_geometries/sim${P}/;
done

echo "Simulations terminated. Log files are in 'Virtual-Angiogenesis/LogFiles/', results in 'Virtual-Angiogenesis/perfused-geometries/sim{x}/'."
