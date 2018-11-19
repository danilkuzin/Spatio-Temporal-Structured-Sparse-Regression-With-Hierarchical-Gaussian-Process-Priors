function [ output ] = logsum( logA, logB )
    
    if(logA > logB)
        output =logA + log(1 + exp(logB-logA));
    else
        output =logB + log(1 + exp(logA-logB));
    end;
     
        
        
