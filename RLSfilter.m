% myRLS.m - recursive least squares algorithm
%
% parameters setting:
%   p     = filter order
%   mu    = step size
%   wn    = fiter weight vector
%   lamda = Exponential weighting factor
%
%Initialization:
%   w0 = 0
%   P = delta^-1 * I
%
%Computation :
%   zn=P*xn
%   k=(1/(lamda+xn'*zn))*(zn)
%   yn=wn'xn
%   en=dn-yn
%   wn+1=wn+k*en
%   Pn+1=(1/lamda)*(Pn-k*zn')
% -------------------------------------------------------
function [en,yn,wn] = RLSfilter(dn,xn,p,lamda)

Length=length(dn);
% xn = xn;
%initialization
I=eye(p);
delta=0.01;
P=1/delta * I;
wn_1=zeros(p,1);
x_1=zeros(p,1);
yn=zeros(Length,1);
en=zeros(Length,1);

%Iteration process
for n = 1:Length
    x_1=[xn(n); x_1(1:p-1)];
    k=(P*x_1)./(lamda+x_1'*P*x_1);
    yn(n)=x_1'*wn_1;
    en(n)=dn(n)-yn(n);
    wn_1=wn_1+k * en(n);
    P=(P-k*x_1'*P)./lamda;
    wn(:,n)=wn_1;
end

end

