function [ output_args ] = Untitled2( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Se piden los datos de cada elemento 
clc;
clear all;
file = input('¿Cual es el nombre del archivo?: ','s');
file = [file,'.xlsx'];
a = xlsread (file);

numerodeelementos = size(a);
%Se organizan las matrices de cada elemento
for i=1:1:numerodeelementos(1,1);
matrizdatos(i,1)= a(i,1);
ditancia (i,1)= a(i,1);

    for j=1:1:numerodeelementos(1,1)  
        if a(i,4)==a(j,6);
            for k=2:1:8
                matrizdatos(i,k)=a(j,k+4);
            end
        elseif a(i,5)==a(j,6);
            for k=9:1:15
                matrizdatos(i,k)=a(j,k-3);
            end
        end
  
    end   
end
%Se organiza las matriz para que quede los elementos de izquierda a derecha
%y de abajo acia arriba.
%for
    
%end

disp ( matrizdatos);
%Calculo de distancias
for i=1:numerodeelementos(1,1)
    ditancia (i,2)= ((matrizdatos(i,5)-matrizdatos(i,12))^2+(matrizdatos(i,6)-matrizdatos(i,13))^2)^0.5;
end







end

