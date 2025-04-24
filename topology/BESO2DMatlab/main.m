clc; clear all; close all;
L=100; %Length
h=50;  %Height
nelx=100;
nely=50;
er=0.01;
volfrac=0.5;
E=1;
nu=0.3;
%f=1;
rmin=3.00;
% rmin=6;
BESO2D(L,h,nelx,nely,volfrac,E,nu,er,rmin)