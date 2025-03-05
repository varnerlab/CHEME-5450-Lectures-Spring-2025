# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# School of Chemical Engineering Purdue University
# W. Lafayette IN 46907 USA
#
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
function Balances(t,x,dxdt_vector,data_dictionary,flow_function)
# ---------------------------------------------------------------------- #
# Balances.jl
# Username: jeffreyvarner
# Type:
# Version: 1.0
# Generation timestamp: 1-19-2016 1:31:08
#
# Arguments:
# t  - current time
# x  - state vector
# dxdt_vector - right hand side vector
# data_dictionary  - Data dictionary instance (holds model parameters)
# ---------------------------------------------------------------------- #

# Correct nagative x's = throws errors in control even if small -
idx = findall(x->(x<0),x);
x[idx] .= 0.0;

# Alias the species vector -
S = x[1];
P = x[2];
X = x[3];
V = x[4];

# Evaluate the reaction kinetics -
(specific_rate_vector) = Kinetics(t,x,data_dictionary);

# Evaluate the flow vector -
(flow_vector) = flow_function(t,x,data_dictionary);

# Evaluate the differential equations -
parameter_array = data_dictionary["KINETIC_PARAMETER_ARRAY"];
YXS = parameter_array[1];
YPS = parameter_array[2];
YXP = parameter_array[3];

# Substrate S - 1
dxdt_vector[1] = flow_vector[1] - ((1/YXS)*specific_rate_vector[1] + (1/YPS)*specific_rate_vector[2])*X;

# Product P - 2
dxdt_vector[2] = flow_vector[2] + specific_rate_vector[2]*X;

# Cellmass X - 3
dxdt_vector[3] = flow_vector[3] + (specific_rate_vector[3] - specific_rate_vector[4])*X;

# Volume V - 4
dxdt_vector[4] = flow_vector[4];

return dxdt_vector;
end
