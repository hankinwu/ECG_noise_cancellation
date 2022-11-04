% myNLMS.m - Normalized least mean squares algorithm
%
% parameters setting:
%   p  = filter order
%   mu = step size
%   wn = fiter weight vector
%   a  = the bias parameter
%   dn = reference signal
%
%Initialization:
%   w0 = 0
%
%Computation :
%   yn=wn'xn
%   k= mu/(a + xn'*xn)
%   en=dn-yn
%   wn+1=wn+k*en*xn
% ------------------------------------------------------------------------
function [en,yn,wn]=NLMSfilter(dn,xn,mu,p,a)

% Initialization
Length=length(dn);
xx=zeros(p,1);
wn_1=zeros(p,1);
yn=zeros(Length,1);
en=zeros(Length,1);

% Iteration process
for n=1:Length
    xx=[xx(2:p);xn(n)];
    yn(n)=wn_1' * xx;
    k=mu/(a + xx'*xx);
    en(n)=dn(n)-yn(n);
    wn_1=wn_1+k*en(n)*xx;
    wn(:,n)=wn_1;
end
end
