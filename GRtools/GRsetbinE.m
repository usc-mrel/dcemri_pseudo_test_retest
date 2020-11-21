function binE = GRsetbinE( tbolus, GRopt, tstamp )

if GRopt.variable_fr && ~isempty(tbolus)    %   Variable frame rate
    
    switch GRopt.app
        case 'perm'
            A = (-1000:1000:1000);
            A = A(A<=-15);
            B = (-1000:5:1000);
            B = B(B>=-15 & B <= +30);
            C = (-1000:10:1000);
            C = C(C>+30 & C <= +80);
            D = (-1000:20:1000);
            D = D(D>+80);
            binE = [A B C D] + tbolus;
        case 'perf'
            A = (-1000:1000:1000);
            A = A(A<=-10);
            B = (-1000:2:1000);
            B = B(B>=-10 & B <= +30);
            C = (-1000:5:1000);
            C = C(C>+30 & C <= +80);
            D = (-1000:20:1000);
            D = D(D>+80);
            binE = [A B C D] + tbolus;
        case 'angio'
            A = (-1000:1000:1000);
            A = A(A<=-2.5);
            B = (-1000:1.25:1000);
            B = B(B>=-2.5 & B <= +15);
            C = (-1000:60:1000);
            C = C(C>=+15 & C <= +80);
            D = (-1000:60:1000);
            D = D(D>+80);
            binE = [A B C D] + tbolus;
        case 'anat'
%             A = -5;
%             A = A(A<=-5);
%             B = (-1000:60:1000);
%             B = B(B>=-5 & B <= +30);
%             C = (-1000:60:1000);
%             C = C(C>+30 & C <= +80);
%             D = (-1000:300:1000);
%             D = D(D>+80);
            %binE = [A B C D] + tbolus;
            %binE = 60:60:1000;
            binE = [-10 30 80 81] + tbolus;

        otherwise
            error('Unrecognized application');
    end
    
    
else    %    Fixed frame rate
    
    if isempty(GRopt.binE)
        FR = 10;
        binE = -FR:FR:1000;
    elseif (length(GRopt.binE)) == 1
        FR = GRopt.binE;
        binE = -FR:FR:1000;
    else
        binE = GRopt.binE;
    end
    
    if ~isempty(tbolus)
        binE = binE + tbolus;
    end
    
end


%   Crop to min/max time
% binE = binE(binE>=0);
binE = binE(binE<max(tstamp));
if binE(end) < tstamp(end)
    binE(end) = tstamp(end);
end

end

