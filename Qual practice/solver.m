clc;
clear;
close all;

lambda = 1;

s = (1:1/1000:15)*lambda;

num = cos(2*pi*s/lambda);

den = sqrt(s);

term = num./den;

exp1 = abs(term-0.577);

figure;plot(s,exp1);