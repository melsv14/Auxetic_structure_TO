clc; clear all; close all;
L=40; %Length
W=10; %Width
h=10;  %Height
nelx=40;
nely=10;
nelz=10;
% BESO3D(40,2,20,20,2,10,0.5,1,0.3,0.01,1.5) % example
% L=20; %Length
% W=20; %Width
% h=20;  %Height
% nelx=2;
% nely=2;
% nelz=2;
er=0.01;
volfrac=0.5;
E=1;
nu=0.3;
%f=1;
rmin=1.25;
BESO3D(L,W,h,nelx,nely,nelz,volfrac,E,nu,er,rmin)