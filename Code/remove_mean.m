function Y = remove_mean(Y, mean_Y)
    Y  = bsxfun(@minus, Y, mean_Y);