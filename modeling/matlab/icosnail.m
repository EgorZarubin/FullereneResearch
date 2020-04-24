clear all;
close all; clc;
k=pi;



for fi=0:0.1:5*pi
 for psi=0:0.1:pi
     r=fi/k;
     x=r*sin(psi)*cos(fi);
     y=r*cos(psi);
     z=r*sin(psi)*sin(fi);
     scatter3(x,y,z,5,'gr'); hold on;
 end
end