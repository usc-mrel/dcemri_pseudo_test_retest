function GRrecon_launch(path,install)

%   Change directory
cd(path);

%   Read pfile header
files = dir('P*.7');
par = read_MR_headers(files(1).name,'all','raw');
bolus_arrival = par.image.user27;
tres = par.image.user28;
flag = par.image.user30;
% [np, nv, ns] = get_key_params(par);
% rp = par.rdb_hdr.rc_xres;
% rv = par.rdb_hdr.rc_yres;



%   Set options
GRopt = GRoptset(flag,'perm');
GRopt.bolus_arrival = bolus_arrival;
if tres > 0
    GRopt.binE = tres:tres:1000;
end
GRopt.DCM_install = str2num(install);
% GRopt.mtrxH = [np 1.1*nv

%   Disable some options
GRopt.variable_fr = 0;

%   Open parallel pool
parpool('NEON',GRopt.parallel,'AttachedFiles',{'startup.m','poolStartup.m','GRrecon_launch.m'});
spmd
   load_TLS;
end


%   Call recon
GRrecon(GRopt,path);

%   Close parallel pool
delete(gcp);
