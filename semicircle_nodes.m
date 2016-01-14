function [ nodes ] = semicircle_nodes( A,B,D,N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

C = (A+B)/2;
CA = A-C; CD = D-C;
r = norm(CA,2);

thita = pi/N;
nodes = zeros(N+1, 3);
nodes(1,:) = A;
nodes(end,:) = B;

nodes(2:end-1, 1) = C(1) + ( CA(1)*cos( thita*(1:N-1)) + CD(1)*sin( thita*(1:N-1)) );
nodes(2:end-1, 2) = C(2) + ( CA(2)*cos( thita*(1:N-1)) + CD(2)*sin( thita*(1:N-1)) );
nodes(2:end-1, 3) = C(3) + ( CA(3)*cos( thita*(1:N-1)) + CD(3)*sin( thita*(1:N-1)) );

end

