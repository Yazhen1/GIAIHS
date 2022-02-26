% Demo
% ========================================================================
% This is the implementation of the algorithm for two-stage PS fusion
% 
% Please refer to the following paper
% 


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




