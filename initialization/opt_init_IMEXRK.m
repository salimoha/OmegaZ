function opt_init_IMEXRK()
global n m ms bnd1 bnd2  Search r MESH_SIZE ConBound RES FUN_MAX CON_MAX DX iter_max
% read the configuration file
conf = readConf('./conf/IMEXRK.conf',0);
n = str2num(conf.NumParam);
ms = str2num(conf.NumConstraints);
m=2*n; % initial starting points
bnd1 = str2num(conf.amin)';
bnd2 = str2num(conf.amax)';
r = str2num(conf.ConstProjc);
Search.method = str2num(conf.method);
Search.constant = str2num(conf.optVar);
iter_max = str2num(conf.MAX_ITER);
MESH_SIZE = str2num(conf.MESH_SIZE);
%  constraints 
ConBound = str2num(conf.constraint_bound);
RES = str2num(conf.delta_tol);
FUN_MAX  = str2num(conf.FUN_MAX);
CON_MAX = str2num(conf.CON_MAX);
DX = str2num(conf.TOL_C); % perturbation for c_i's
% keyboard
end

