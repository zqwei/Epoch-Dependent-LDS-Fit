function Y = submat (Y, rowmin, rowmax, colmin, colmax)

if isempty(rowmin); rowmin = 1; end
if isempty(rowmax); rowmax = size(Y,1); end
if isempty(colmin); colmin = 1; end
if isempty(colmax); colmax = size(Y,2); end

Y = Y(rowmin:rowmax, colmin:colmax);