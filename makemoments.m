function [mom] = makemoments(par,sim)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

disp("Calculate moments...")
disp("Not implemented yet")

mom = struct();

ave_prof  = mean(sim.prof(:));   % mean profitability
ave_irate = mean(sim.Irate(:)); % mean investment
std_irate = std(sim.Irate(:));  % volatility of investment

mom.ave_prof  = ave_prof;
mom.ave_irate = ave_irate;
mom.std_irate = std_irate;


end

