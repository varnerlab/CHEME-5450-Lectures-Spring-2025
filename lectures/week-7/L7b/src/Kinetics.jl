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
function Kinetics(t,x,data_dictionary)

  # Alias the species vector -
  S = x[1];
  P = x[2];
  X = x[3];

  # We have four rates, qS, qP, mu and death -
  specific_rate_vector = zeros(4);

  # Get parameters from the data file -
  parameter_array = data_dictionary["KINETIC_PARAMETER_ARRAY"];
  YXS = parameter_array[1];
  YPS = parameter_array[2];
  YXP = parameter_array[3];
  mugmax = parameter_array[4];
  KGS = parameter_array[5];
  kd = parameter_array[6];
  ms = parameter_array[7];
  mp = parameter_array[8];

  # calculate the specific growth rate -
  specific_grow_rate = mugmax*(S)/(KGS + S);

  # specific_rate_vector[1] = qS
  if (specific_grow_rate<0.1*mugmax)
    specific_rate_vector[1] = specific_grow_rate;
  else
    specific_rate_vector[1] = specific_grow_rate - YXS*ms;
  end

  # specific_rate_vector[2] = qP  
  if (specific_grow_rate<0.1*mugmax)
    specific_rate_vector[2] = (1/YXP)*specific_grow_rate + mp;
  else
    specific_rate_vector[2] = (1/YXP)*specific_grow_rate;
  end

  # specific_rate_vector[3] = mu
  specific_rate_vector[3] = specific_grow_rate;

  # specific_rate_vector[4] = specific death rate
  specific_rate_vector[4] = kd;

  # return the rate vector -
  return specific_rate_vector;
end
