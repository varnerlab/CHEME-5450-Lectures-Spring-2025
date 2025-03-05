include("Balances.jl")
using Sundials;

# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# School of Chemical Engineering Purdue University
# W. Lafayette IN 46907 USA

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
function SolveBalances(TSTART,TSTOP,Ts,data_dictionary,flow_function)
# ----------------------------------------------------------------------------------- #
# SolveBalances.jl
# SolveBalances: Solves model equations from TSTART to TSTOP given parameters in data_dictionary.
# Username: jeffreyvarner
# Type:
# Version: 1.0
# Generation timestamp: 1-19-2016 1:33:16
#
# Input arguments:
# TSTART  - Time start
# TSTOP  - Time stop
# Ts - Time step
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# TSIM - Simulation time vector
# X - Simulation state array (NTIME x NSPECIES)
# ----------------------------------------------------------------------------------- #

# Get required stuff from DataFile struct -
TSIM = collect(TSTART:Ts:TSTOP);
initial_condition_vector = data_dictionary["INITIAL_CONDITION_ARRAY"];

# Call the ODE solver -
fbalances(t,y,ydot) = Balances(t,y,ydot,data_dictionary,flow_function);
X = Sundials.cvode(fbalances,initial_condition_vector,TSIM,reltol=1e-8,abstol=1e-10);

return (TSIM,X);
end
