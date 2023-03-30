function [ X ] = AnyRand( PDF, a, b )
%function [ X ] = AnyRand( PDF, a, b )
%   Generates ONE random number of any PDF froma a to b. To get a sample
%   you must loop over this function the desired number of elements.

pdfmax = 1;
integers = 0;
maxite = 10000;

for i = 1:1:maxite
    if integers == 1
        % not yet
    else
        rand_x = (b - a) * rand(1) + a;

    rand_y = pdfmax * rand(1);
    calc_y = PDF(rand_x);

    if (rand_y <= calc_y)
        X =  rand_x;
        return
    end
    end
end

if i == maxite
    error(['Could not find a matching random number within pdf in',...
           string(maxite) + 'iterations.'])
end


end

