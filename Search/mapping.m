function [xp,xq,isxqp] = mapping(SolC)
% This function is for mapping the c_i in IMEXRK project
% x: quantized SolC
% x1: x perturb
conf = readConf('./conf/IMEXRK.conf',0);
N = str2num(conf.MESH_SIZE);
m = length(SolC);
x=round(SolC'*N);
[xa,ind]=sort(x);
Ain=[eye(m) zeros(m); -eye(m) zeros(m); eye(m) -eye(m); -eye(m) -eye(m); 1 -1 0 zeros(1,m); 0 1 -1 zeros(1,m)];
bin=[N*ones(m,1); zeros(m,1); xa; -xa; -1; -1];
f=[zeros(m,1) ;ones(m,1)]; 
options = optimoptions('linprog','Algorithm','dual-simplex','display','none');
x1 = linprog(f,Ain,bin,[],[],[],[],[],options);
if sum(x1(m+1:end)) < 1e-10
    isxqp = 1;
else
    isxqp = 0;
end
x1=x1(1:m);
x1(ind)=x1;
xp=x1./N;
xq=x./N;
end