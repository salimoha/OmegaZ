function [sigma,mu,transtime]= initial_cal(x,index)
% Calculate an initial calculation at a new point

% make the folder: 00001, 

% generate the gometry: python file, put the geometry (.ply) in the folder  

% Calculate the text files and initial time series: ?????

% run NFA:

% Read the time series file
[t,zf]=Read_log_file('Zforces.log');
% Calculate the statitical analysis
Initial_mod(1)=0; Initial_mod(2)=0;
[mu,sigma,Initial_mod] = alpha_statitical_analyze(zf,Initial_mod);
transtime=Initial_mod(2);
end