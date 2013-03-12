function [y, maxx] = lowexp(x)
%Takes the exp of some logs, but adds a const factor to prevent Infs and
%NaNs.

%Richard Mann (2010)
maxx = max(x(:));
y = exp(x-maxx);



