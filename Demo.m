% Demo
% ========================================================================
% This is the implementation of the algorithm for two-stage PS fusion
% 
% Please refer to the following paper
% Wang, Y.; Liu, G.; Zhang, R.; Liu, J. A Two-Stage Pansharpening Method for the Fusion of Remote-Sensing Images. Remote Sens. 2022, 14, 1121. https://doi.org/10.3390/rs14051121 


clc; 
clear
close all;

%==========================================================================
disp('==============================================================');
disp('In order to illustrate the problem');
disp('the scene we selected contains large areas of vegetation!');
fprintf('...\n');
load('Gaofen.mat')
Pan=double(imgPAN);
Mul11 = double(imgMS);Mul = imresize(Mul11,size(Pan),'bicubic');
F4= GIAIHS(Mul,Pan);
disp('Done!')
disp('==============================================================');




