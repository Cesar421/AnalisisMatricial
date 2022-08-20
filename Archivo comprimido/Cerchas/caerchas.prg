function [ output_args ] = parr( input_args )
 

%CALCULO DE REACCINES cerchass
%este programa analiza parrillas con acciones exclusivamente sobre el eje z respecto a un plano
%coordenado global.


%clc;
 
disp ('Recuerde manejar siempre las mismas unidades en el ingreso de datos');
disp ('Ademas recuerde ingresar las coordenadas de los nodos SIEMPRE de izquierda a derecha del menor a mayor punto');
disp ('cuando ingrese los grados de libertad se empiesa siempre a enumerar de acuerdo a los desplazamientos desconocidos')
disp ('  ');
disp ('  ');
 
disp ('                                 INICIO        ');
disp ('Programa que calcula reacciones en una parrilla con acciones solo en el eje z ');
 
barrasnum =  input ('digite el número de barras    ');
nodos = input ('digite el número de nodos    ');
dc = input ('digite el número de desplazamientos conocidos    ');
GIC = 2*nodos-dc;
disp (' ');
disp ('grado de indeterminación cinematico    ');
fprintf ('GIC = %i.  ',GIC); 
disp ('   ');