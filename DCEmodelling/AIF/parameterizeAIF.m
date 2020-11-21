function [ param ] = parameterizeAIF( AIF, opt )
%function [ param ] = parameterizeAIF( AIF, opt )
%   fits parameterized model to AIFs
%
%   this is just a wrapper for 'fitAIF' to generate a good initial guess
%
%   Input:
%       AIF                     concentration time curve
%       opt                     various options
%       .fid                    file handle to output
%       .time                   time vector
%       .AIF.paramterization    currenlty only '2Gammas+Sigmoid'
%       .AIF.initialParameters  container of initial parameters
%                                 for exponential sigmoid (2):
%                                   .a      amplitudes [2,1]
%                                   .delta 	time delays [2,1]
%                                   .alpha 	upslope parameter [2,1]
%                                   .tau   	downslope parameter [2,1]
%                                   .beta   sigmoid parameter [1,1]
%   Output:
%       same as initialParam
%
% Yannick 2018
%

switch opt.AIF.paramterization
    
    case '2Gammas+Sigmoid'
        if isempty(opt.AIF.initialParameters)
            timeScale = opt.time(2) - opt.time(1);
            
            [~, Imax] = max(AIF);
            
            localMin = (AIF(1:end-2) > AIF(2:end-1)) & (AIF(3:end) > AIF(2:end-1));
            localMin = [false localMin false];
            localMin = find( localMin & (1:length(AIF) > Imax) & ((opt.time - opt.time(Imax) < 60)) );
            
            localMax = (AIF(1:end-2) < AIF(2:end-1)) & (AIF(3:end) < AIF(2:end-1));
            localMax = [false localMax false];
            localMax = find( localMax & (1:length(AIF) > Imax) & ((opt.time - opt.time(Imax) < 40)) );
            
            if length(localMin) < 1
                localMin = Imax + 2;
            end
            
            if length(localMax) < 1
                localMax = localMin(1) + 1;
            end
            
            if length(localMin) < 2
                localMin(end+1) = localMax(1) + 1;
            end
            
%             if length(localMax) < 1
%                 recircBump = false;
%             else
%                 % de-trend the tail of AIF to estimate the noise power and
%                 % figure out if the recirculation bump is real
%                 AIFtailStart = localMin(2);
%                 AIFtail = AIF(AIFtailStart:end);
%                 p = polyfit(1:length(AIFtail),AIFtail,2);
%                 AIFtail = AIFtail - polyval(p, 1:length(AIFtail));
%                 estBump = AIF(localMax(1)) - polyval(p, localMax(1));
%                 noisePower = std(AIFtail);
%                 
%                 % So, yay or nay on the recirculation bump?
%                 recircBump = estBump > 2*noisePower;
%             end
            

            
            % it is important for the initial guess to have the peak
            % correct!
%             if recircBump
            initParam = struct();
            
            initParam(1).a     = [max(AIF(1:localMin(1))) max(AIF(localMin(1):localMin(2)+1)) max(AIF(localMin(2):end))];
            initParam(1).alpha = [3  3 0.015 / timeScale];
            initParam(1).tau   = [0.9 0.9 5] * timeScale;
            initParam(1).delta = [max(0,Imax-3) localMax(1)-2 localMin(2)-4] * timeScale + opt.time(1);
            initParam(1).beta  = 1.3 / timeScale;
            num_param(1) = 13;
%             else
            initParam(2).a     = [max(AIF(1:localMin(1))) max(AIF(localMin(1):end))];
            initParam(2).alpha = [3  0.015 / timeScale];
            initParam(2).tau   = [0.9 5] * timeScale;
            initParam(2).delta = [max(0,Imax-3) max(Imax, localMin(1)-6)] * timeScale + opt.time(1);
            initParam(2).beta  = 1.3 / timeScale;
            num_param(2) = 9;
%             end
        else
            initParam = opt.AIF.initialParameters;
        end
        
        SC = inf(1,length(initParam));
        for i=1:length(initParam)
            [ param(i) ] = fitAIF( AIF(:), opt.time, 2, initParam(i) );
            fittedAIF = AIFgammavariates( opt.time, param(i) );
            SSE = norm(AIF(:) - fittedAIF(:));
            SC(i) = AIC_gaussian(num_param(i) + 1, length(AIF), SSE, 1);
        end
%         
%         close all
%         plot( opt.time, AIF), hold on
%         for i=1:length(initParam)
%             initAIF = AIFgammavariates( opt.time, initParam(i) );
%             fittedAIF = AIFgammavariates( opt.time, param(i) );
% 
%             plot( opt.time, initAIF, '-', 'linewidth', 1.5 )
%             plot( opt.time, fittedAIF, '--', 'linewidth', 1.5 )
%         end
        
        [~, I] = min(SC);
        param = param(I);

    otherwise
        err_msg = 'Unkown AIF model to fit to!';
        fprintf(opt.fid, '\t\t\t ERROR: %s \n', err_msg);
        error(err_msg)
end

end