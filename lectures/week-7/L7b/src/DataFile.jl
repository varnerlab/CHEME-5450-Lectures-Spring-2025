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
function DataFile(TSTART,TSTOP,Ts)
# ----------------------------------------------------------------------------------- #
# DataFile.jl
# DataFile: Stores model parameters as key - value pairs in a Julia Dict()
# Username: jeffreyvarner
# Type:
# Version: 1.0
# Generation timestamp: 1-19-2016 1:34:33
#
# Input arguments:
# TSTART  - Time start
# TSTOP  - Time stop
# Ts - Time step
#
# Return arguments:
# data_dictionary  - Data dictionary instance (holds model parameters)
# ----------------------------------------------------------------------------------- #

# Initial conditions in the reactor -
initial_condition_array = zeros(4);
initial_condition_array[1] = 10.0;  # 1 S
initial_condition_array[2] = 0.0;  # 2 P
initial_condition_array[3] = 0.1;  # 3 X
initial_condition_array[4] = 1.0;  # 4 V


# Kinetic parameters -
kinetic_parameter_array = zeros(8);
kinetic_parameter_array[1] = 0.5;   # 1 YXS
kinetic_parameter_array[2] = 0.6;   # 2 YPS
kinetic_parameter_array[3] = 0.5;   # 3 YXP
kinetic_parameter_array[4] = 0.75;   # 4 mugmax
kinetic_parameter_array[5] = 1.5;   # 5 KGS
kinetic_parameter_array[6] = 0.05;   # 6 kd
kinetic_parameter_array[7] = 0.08;   # 7 maintenance
kinetic_parameter_array[8] = 0.0;   # 8 non-specific production

# Feed parameters -
feed_parameter_array = zeros(3);
feed_parameter_array[1] = 50.0;   # 1 SIN
feed_parameter_array[2] = 0.0;    # 2 PIN
feed_parameter_array[3] = 0.0;    # 3 XIN

# For chemostats, we set the diltuin rate -
dilution_rate = 0;

# For fed-batch culutes we set parameters which appear in the
# input flow model -
flow_model_parameter_array = zeros(3);
flow_model_parameter_array[1] = 0.2;
flow_model_parameter_array[2] = 0.005;
flow_model_parameter_array[3] = 5.0;

# ---------------------------- DO NOT EDIT BELOW THIS LINE -------------------------- #
data_dictionary = Dict();
data_dictionary["FEED_PARAMETER_ARRAY"] = feed_parameter_array;
data_dictionary["KINETIC_PARAMETER_ARRAY"]= kinetic_parameter_array;
data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;
data_dictionary["DILUTION_RATE"] = dilution_rate;
data_dictionary["FLOW_MODEL_PARAMETER_ARRAY"] = flow_model_parameter_array;
# ----------------------------------------------------------------------------------- #
return data_dictionary;
end
