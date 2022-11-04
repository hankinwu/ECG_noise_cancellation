% LMSfilter.m - least mean squares algorithm
%
% parameters setting:
%   p  = filter order
%   mu = step size
%   wn  = fiter weight vector
%
%Initialization:
%   w0 = 0
%
%Computation :
%   yn=wn'xn
%   en=dn-yn
%   wn+1=wn+mu*en*xn
%
% -------------------------------------------------------
function [en,yn,wn] = LMSfilter(dn,xn,mu,p)
Length=length(dn);

% Initialization
xx=zeros(p,1);
wn_1=zeros(p,1);
yn=zeros(Length,1);
en=zeros(Length,1);

% Iteration process
for n=1:Length
    xx=[xx(2:p);xn(n)];
    yn(n)=wn_1'*xx;
    en(n)=dn(n)-yn(n);
    wn_1=wn_1+mu*en(n)*xx;
    wn(:,n)=wn_1;
end

end

