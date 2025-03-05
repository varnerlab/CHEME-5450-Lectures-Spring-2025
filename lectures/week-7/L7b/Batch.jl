# Script to solve the simple model -
include("Include.jl");

# ================================================================== #
# Phase 0 = Initialize
# ================================================================== #
data_dictionary = DataFile(0,0,0);
flow_function = BatchFlow;

# ================================================================== #
# Phase 1 = run for 12 hours
# ================================================================== #

# Setup the timescale for phase 1 -
TSTART_P1 = 0;
TSTOP_P1 = 12.0;
Ts_P1 = 0.01;

# Solve the balance equations -
(TP1,XP1) = SolveBalances(TSTART_P1,TSTOP_P1,Ts_P1,data_dictionary,flow_function);
