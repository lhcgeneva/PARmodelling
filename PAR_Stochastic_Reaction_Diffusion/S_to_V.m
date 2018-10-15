function [StoV, S, V, outpt1, outpt2] = S_to_V( mode, describers )
%S_TO_V Calculates surface to volume ratio...
%   ...for a prolate spheroidal geometry, outpt1 and outpt2 depend on the
%   mode chosen see condition at the very bottom.

if strcmp(mode, 'Circumference')
    circ = describers{1};
    aspRatio = describers{2}; % Short over long axis
    e = (1 - aspRatio.^2);  %No square because of definition of parameter in 
                            %ellipke() without square!
    [ ~, EllInt2nd] = ellipke(e);
    longAxis = circ./(4*EllInt2nd);
    shortAxis  = longAxis .* aspRatio;
elseif strcmp(mode, 'Axes')
    shortAxis = describers{1};
    longAxis  = describers{2};
elseif strcmp(mode, 'Volume')
    vol = describers{1};
    aspRatio = describers{2};
    shortAxis = (3*aspRatio.*vol/(4*3.1415)).^(1/3);
    longAxis = shortAxis./aspRatio;
else
    disp('Error in S_to_V: Mode not supported');
    return
end


e = sqrt(1 - shortAxis.^2./longAxis.^2);
S = 2 * 3.1415 * shortAxis.^2 .* (1 + longAxis./(shortAxis .* e) .* asin(e));
V = 4/3 * 3.1415 * shortAxis.^2 .* longAxis;
StoV = S./V;

if strcmp(mode, 'Circumference')
    outpt1 = shortAxis;
    outpt2 = longAxis;
elseif strcmp(mode, 'Axes')
    outpt1 = shortAxis./longAxis;
    [~, e] = ellipke((1-shortAxis.^2./longAxis.^2));
    outpt2 = 4. * longAxis .* e; %Circumference
elseif strcmp(mode, 'Volume')
    outpt1 = shortAxis;
    [~, e] = ellipke((1-shortAxis.^2./longAxis.^2));
    outpt2 = 4. * longAxis .* e; %Circumference
else
    disp('Haha, you can''t even get here! :P');
end

end

