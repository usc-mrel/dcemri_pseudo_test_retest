%   Launch parallel pool
parpool(12);

%   Run anatomic recon
t = tic;
GRopt = GRoptset(2,'anat');
GRrecon(GRopt);

%   Run angiography recon
GRopt = GRoptset(26,'angio');
GRrecon(GRopt);

%   Run perfusion recon
GRopt = GRoptset(8,'perf');
GRrecon(GRopt);

%   Run permeability recon
GRopt = GRoptset(8,'perm');
GRrecon(GRopt);

t = toc(t);
fprintf('Time = %g (hours)\n',t/3600);

%   Close parallel pool
delete(gcp);

