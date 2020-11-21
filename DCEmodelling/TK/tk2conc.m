function conc = tk2conc(TK, AIF, opt)
%function conc = tk2conc(TK, AIF, opt)
%
%   compute concentration time curves from TK parameters and AIF
%
% Yannick Bliesener 2019


N = prod(opt.size(1:3));
Nt = opt.size(4);

% computed shifted AIFs
Nshifts = 1;
if ~isempty(opt.AIF.shifts)
    Nshifts = length(opt.AIF.shifts);
    shiftedAIFs = zeros(Nt, Nshifts);
    for s=1:Nshifts
        dtime = opt.time(2)-opt.time(1);
        shiftedAIFs(:,s) = shiftAIF(AIF, opt.time, opt.AIF.shifts(s)*dtime);
    end
else
   shiftedAIFs = AIF; 
end

% choose last model (assuming it is the most complex model) as overall
% model
tkmodel = opt.TK.models{end};

% compute contration time curves
conc = zeros([N Nt]);
for s=1:Nshifts
    n=find(TK.shift == s);
    switch tkmodel
        case 'null'
            % pass
        case 'vpmodel'
            intAIF = integrateAIF(opt.time, shiftedAIFs(:,s), opt.AIF.integration);                
            conc(n,:) = model_patlak(length(n), Nt, TK.vp(n), zeros(length(n),1), shiftedAIFs(:,s), intAIF);
        case 'patlak'
            intAIF = integrateAIF(opt.time, shiftedAIFs(:,s), opt.AIF.integration);                
            conc(n,:) = model_patlak(length(n), Nt, TK.vp(n), TK.kt(n), shiftedAIFs(:,s), intAIF);
        case 'etk'
            conc(n,:) = model_standard(length(n), Nt, TK.vp(n), TK.kt(n), TK.kep(n), shiftedAIFs(:,s), opt.time, opt.AIF.integration);
        case 'tcxm'
            conc(n,:) = model_2cxm(length(n), Nt, TK.vp(n), TK.ve(n), TK.kt(n), TK.fp(n), shiftedAIFs(:,s), opt.time, opt.AIF.integration);
        otherwise
            error('Unknown model to fit concentration time curves')
    end
end

conc = reshape(conc, opt.size(1:4));

end

