
%CALCULO DE REACCINES cerchass
%clc;
 
disp ('Recuerde manejar siempre las mismas unidades en el ingreso de datos');
disp ('Ademas recuerde ingresar las coordenadas de los nodos SIEMPRE de izquierda a derecha del menor a mayor punto');
disp ('cuando ingrese los grados de libertad se empiesa siempre a enumerar de acuerdo a los desplazamientos desconocidos')
disp ('  ');
disp ('  ');
 
disp ('                                 INICIO        ');
disp ('Programa que calcula reacciones en una cercha');
c=  input ('digite el número de rotulas   '); 
barrasnum =  input ('digite el número de barras    ');
nodos = input ('digite el número de nodos    ');
dc = input ('digite el número de desplazamientos conocidos  ');
GIE=barrasnum+dc-2*nodos-c;
GIC = 2*nodos-dc;
disp (' ');
disp ('grado de indeterminación cinematico    ');
fprintf ('GIC = %i.  ',GIC); 
disp ('grado de indeterminación estatico   ');
GIE
disp ('   ');

for cont=1:barrasnum ;
    
    %--------------------------------------------------------- 
    %podria crearse un menu que permita que el usuario diga si todos los elementos de la estructura
    %son iguales de esta manera se Ahorra tiempo a la hora de ingresar los datos.
    
  % se almacenan las propiedades de cada barra  tanto geometricas como fisicas, ademas tambien se solicita las cargas que hay en este elemento 
 
     disp (' ');
     fprintf ('PROPIEDADES GEOMETRICAS DEL ELEMENTO NUMERO %i.    ',cont );
     disp (' ');
     disp (' ');
    
     disp ('1)  seccion transversal circúlar');
     disp ('2)  seccion transversal rectangular');
     disp ('                 ');
     a = input ('seleccione un opción    ');
             
     %sección circular
     if a == 1;
        disp(' ');
        barras{cont}.radio = input('digite el radio de la seccion transversal del elemento    '); 
        barras{cont}.I= (3.1416*(barras{cont}.radio)^4)/4;
        barras{cont}.J = (3.1416*(barras{cont}.radio)^4)/2;
        barras{cont}.area=pi()*(barras{cont}.radio)^2
        barras{cont}.L=input('digite la longitud del elemento')
        disp (' ');
        
     %seccion rectangular
     
     elseif a == 2
         
        disp (' ');
        barras{cont}.base = input('digite la base de la seccion transversal del elemento    ');
        barras{cont}.altura = input('digite  la altura de la seccion transversal del elemento    ');
        barras{cont}.I=((barras{cont}.base*(barras{cont}.altura)^3)/12);
        barras{cont}.L=input('digite la longitud del elemento')
        % este condicional funciona para la formula de J recordar que en
        % esta formula t es el valor mas pequeño y b es el mayor esto
        % cambia el resultado de la formula J
        
        if barras{cont}.base >= barras{cont}.altura;
            
        b = barras{cont}.base;
        t = barras{cont}.altura;
        barras{cont}.J=(1/3-0.21*((t/b)*(1-(1/12)*(t/b)^4)))*b*(t)^3;
        barras{cont}.area=b*t
        elseif barras{cont}.altura > barras{cont}.base;
            
        t=barras{cont}.base;
        b=barras{cont}.altura;
        barras{cont}.J=((1/3)-0.21*(t/b)*(1-(1/12)*(t/b)^4))*b*(t)^3;
        barras{cont}.area=b*t    
        end
    
        %disp (barras{cont}.J);% por si se necesita 
        %verificador si daba al J disp (barras{cont}.J) SI DA;
        
        disp (' ');
                
     end
     
     
     fprintf ('PROPIEDADES DEL MATERIAL DEL ELEMENTO NUMERO %i.   ',cont );
     disp (' ');
     disp (' ');
     
     barras{cont}.E = input ('digite el modulo de elasticidad de este elemento (E)    '); 
     disp (' '); 
     angulo=input ('digite el angulo formado el eje local al global en sentido antihorario')
     rad=(angulo*pi()/180)
     A=barras{cont}.area
     barras{cont}.cos=cos(rad);
     barras{cont}.sen=sin(rad);
     C=barras{cont}.cos
     S=barras{cont}.sen
     barras{cont}.Krigidez=(A*barras{cont}.E/barras{cont}.L)*[C^2 C*S -C^2 -C*S;
                                                              C*S S^2 -C*S -S^2;
                                                            -C^2 -C*S C^2 C*S;
                                                            -C*S -S^2 C*S S^2]
    K=barras{cont}.Krigidez
    sizeK=size(K)
    
      % el siguiente for le permite al usuario la ubicación de los grados de
    % libertad de la estructura de izquierda a derecha de menor a mayor
    
    disp ('');
    fprintf ('ASIGNACION DE LOS GRADOS DE LIBERTAD DEL ELEMENTO NUMERO %i:',cont);
    disp (' ');
    disp (' ');
    disp ('por favor digite el valor del grado de libertad de cada nudo del elemento en el siguiente orden:');
    disp ('Xi, Yi,Xj,Yj ');
    disp (' ');
    
    for x = 1 : 4;
        
        barras{cont}.gradosdelibertad(x)=input('');
    end
end
    
    Fext=zeros(nodos*2-dc,1)
     disp (' ');
    numf = input ('¿ Cuantas fuerzas actuan en los nodos de la estructura ?.    ');
    Fext=zeros(nodos*2-dc,1)
 
for x = 1 :numf
    
    disp (' ');
    fprintf ('Indique el grado de libertad donde se ubica la fuerza numero %i.',x)
    disp (' ');
    disp (' ');
    psfuervec = input ('')
    
    disp (' ');
    fprintf ('Ingrese el valor de la fuerza que actúa en el grado de libertad numero %i.    ',psfuervec)
    disp (' ');
    disp (' ');
    funu = input (' ')
    
    Fext(psfuervec)= funu;
    
end
% vector desplazamientos conocidos
 
disp (' ');
dc0=input('¿ Existe algún desplazamiento conocido diferente de cero ? [Y][N]:    ','s');
 
if dc0==('Y'|'y')
    
    disp (' ');
    DB=zeros(dc,1);
    numdcd0=input('¿Cuántos ?');
    
    for x=1:numdcd0
        
         posivec=input('Por favor indique el grado de libertad del desplazamiento conocido diferente de cero:    ');
         desplaz=input('Por favor indique el valor del desplazamiento:    ');
         DB(posivec-(nodos*2-dc))=desplaz;
    end
else
    DB=zeros(dc,1);
end
KGENERAL=zeros(nodos*2,nodos*2);
 
for cont = 1 : barrasnum;
    
    barras{cont}.NUMERODDEGDL =size(barras{cont}.Krigidez);
    
    for i=1 : barras{cont}.NUMERODDEGDL;
        
        ifil = barras{cont}.gradosdelibertad(i); 
        
        for j=1 : barras{cont}.NUMERODDEGDL;
        
            icol = barras{cont}.gradosdelibertad(j);
            
            KGENERAL(ifil , icol) = KGENERAL(ifil, icol) + barras{cont}.Krigidez(i,j) ;
            
        end
        
    end
    disp('  ')
   KGENERAL 

end
 %%se copia la parte RAA de la matriz de rigidezes.
 
RAA =zeros (nodos*2-dc,nodos*2-dc);
 
for i = 1 : nodos*2-dc
    
    for j =1 : nodos*2-dc
        
        RAA (i,j) = KGENERAL (i,j);
    end
    
end
RAA 
%se copia la parte RAB  de la matriz de rigidezes.
 
RAB =zeros(nodos*2-dc,dc);
 
x = 1; % contadores 
y = 1;
 
for i = 1:nodos*2-dc
    
    for j = (nodos*2-dc + 1):nodos*2
        
        RAB(x,y) = KGENERAL (i,j);  
        
    end
    
    x=x+1;
    y=1;
end
RAB 
% se copia la parte de RBA de la matriz de rigidezes.
 
RBA = zeros (dc,nodos*2-dc);
 
x = 1;
y = 1;
 
for i = (nodos*2-dc + 1):nodos*2
    
    for j = 1 : (nodos*2-dc)
        
        RBA (x,y) = KGENERAL(i,j);
        y = y+1;
    end
    y = 1;
    x= x+1;
    
end

RBA
% se copia la parte de RBB de la matriz de rigidezes.
 
RBB = zeros(dc,dc);
 
x = 1;
y = 1;
 
for i=(nodos*2-dc+1):nodos*2
    
    for j = (nodos*2-dc+1):nodos*2
        
    RBB (x,y) = KGENERAL(i,j);     
        y = y + 1;
        
    end
    
    x = x + 1;
    y =1;
    
end
RBB
% se plantean y se resuelven las ecuaciones para encontrar los
% desplazamientos desconocidos y las fuerzas desconocidas

deltasdesconocidos = (inv(RAA))*(Fext-RAB*DB)
disp('Reacciones en los apoyos  ') 
ABN = RBA*deltasdesconocidos  + RBB*DB
Internas=0;
Internas=input('desea conocer las fuerzas internas digite Y/y   ')

if Internas==('y'|'Y')
    disp('digite el grado de libertad de los nudos libres')
    

else
    disp('......fin....')
end 

%%-------------------------------------------------------------------
