%% CONTROL THEORY INTRO - HW1%%
clc;
clear all;
close all;

% Q2: plotting step response and zero-poles map -

A_21 = [1 1];
B_21 = [1 -1 6];
plot_zeros_poles(A_21,B_21, 2.1);

A_22 = [2 1 -1];
B_22 = [3 -2 6];
plot_zeros_poles(A_22,B_22, 2.2);

A_23 = [1 6];
B_23 = [1 3 2];
plot_zeros_poles(A_23,B_23, 2.3);

A_24 = [1];
B_24 = [1 4 4];
plot_zeros_poles(A_24,B_24, 2.4);


% Q3: plotting step response and zero-poles map -

A_31 = [1];
B_31 = [1 1 3 1 1];
plot_zeros_poles(A_31,B_31, 3.3);

% Q4: plotting step response and zero-poles map -

A_41 = [1 0]; %% NEED TO CHANGE ---------------------------------------------
B_41 = [1 1 3 1 1];
plot_zeros_poles(A_41,B_41, 4.3);

% Q5: plotting step response for different times -
T = [2 5 10];
A_i = [1];
plot_5(T, A_i);

% Q5: plotting zero-poles map -
A_61 = [5000];
B_61 = [1 20 1000 5000];
plot_6(A_61, B_61);

% FUNCTIONS ---------------------------------------------------------------

% Q1: function -
function [poles, zeros] = pzplot2(a, b)
    % input: coefficient vector of a transfer function
    % output: plots the poles-zeros map (and highlights the ROC)
    % a: zeros, b: poles
    
    transfer = tf(a,b); % creates transfer function
    zeros = roots(a); % matlab function (recommended)
    poles = roots(b);
    xmax = max([real(zeros);real(poles)]);
    xmin = min([real(zeros);real(poles)]);
    ymax = max([imag(zeros);imag(poles)]);
    ymin = min([imag(zeros);imag(poles)]);
    limits = [xmin-1, xmax+1, ymin-1, ymax+1]; % making sure the axis are at least longer than 1

    % plot figure - 
    pzplot(transfer); % matlab function (recommended)
    axis(limits);
    Real_ROC = [max(real(poles)),xmax+1,xmax+1,max(real(poles))];
    Imaginary_ROC = [ymax+1,ymax+1,ymin-1,ymin-1];
    patch(Real_ROC,Imaginary_ROC, 'blue');
    alpha(0.2);
    xlabel('Real (\sigma)')
    ylabel('Imaginary (j\omega)')
    xline(0);
    yline(0);
end

% (function for Q2 - )
function plot_zeros_poles(a, b, num)
    transfer = tf(a,b); % matlab function (recommended)
    figure;
    subplot(2,1,1);
    stepplot(transfer); % matlab function (recommended) - does both step and plot
    title(['Step Response - Question number: ', num2str(num)])
    hold on;
    subplot(2,1,2);
    pzplot2(a,b);
    title(['Zero-Poles map - Question number: ', num2str(num)])
    hold off;
end

function plot_5(times_list, A)
    figure
    colours = ['g', 'b', 'r'];
    for index = 1:length(times_list)
        temp_num = 5 + (index *0.1);
        B_i = [times_list(index) 1];
        transfer = tf(A, B_i);
        stepplot(transfer,colours(index));
        hold on
    end
    title(['Step Response - Question number: ', num2str(5.1)])
    legend('T=2','T=5','T=10')
    hold off;

end

% -------------------------------------------------------------------------
% 6.2
function plot_6(A_list, B_list)
    C_list = [B_list, 0];
    [r,p,k] = residue(A_list,C_list);
    disp(r);
    disp(p);
    disp(k);
    figure;
    pzplot2(A_list,B_list);
    title('Zero-Poles map - Question number: 6.1')

    % calc y(t) 
    dt = 1/1000;
    t = 0:dt:2;
    y = zeros(size(t));
    y = r(1)*exp(p(1)*t)+r(2)*exp(p(2)*t)+r(3)*exp(p(3)*t)+r(4)*exp(p(4)*t);
    figure
    plot(t,y)
    xlabel('t (sec)')
    ylabel('y(t)')
    title('y(t) - Question number: 6.4')
end

% -------------------------------------------------------------------------
