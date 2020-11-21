
hans = 0;
multi = 1;


if ~hans    
    if multi
        
        mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=10 " -lpthread fwt1dMEX.c
        mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=10 " -lpthread iwt1dMEX.c              
    else
        
        mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=1 " -lpthread fwt1dMEX.c
        mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=1 " -lpthread iwt1dMEX.c        
    end
else
    if multi
        
        mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=12 " -lpthread fwt1dMEX.c
        mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=12 " -lpthread iwt1dMEX.c              
    else
        
        mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=1 " -lpthread fwt1dMEX.c
        mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=1 " -lpthread iwt1dMEX.c        
    end
end


clear hans multi;





