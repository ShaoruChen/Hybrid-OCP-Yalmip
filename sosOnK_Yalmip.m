function [F,c,S] = sosOnK_Yalmip(F, p, x, h, d, c)
%SOSONK Summary of this function goes here
%   Detailed explanation goes here

    m = size( h, 1 ); % total number of h's
    S = sdpvar(m,1);
%     bases = monomials( x, 0:d );
    for i = 1:m
%         [stemp,ctemp,~] = polynomial(x,d,d);
        [stemp,ctemp,~] = polynomial(x,d,d);

        S(i) = stemp;
        c = [c; ctemp];
        F = [F, sos(S(i))]; % ensure that the polynomial is SOS
    end
    [ F ] = [F, sos( p - S'*h )]; % ensure that p is SOS on K 


end
