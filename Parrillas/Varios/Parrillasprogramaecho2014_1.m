function [ output_args ] = parr( input_args )
%CALCULO DE REACCINES PARRILLAS
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
GIC = 3*nodos-dc;
disp (' ');
disp ('grado de indeterminación cinematico    ');
fprintf ('GIC = %i.  ',GIC); 
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
        disp (' ');       
     %seccion rectangular    
     elseif a == 2        
        disp (' ');
        barras{cont}.base = input('digite la base de la seccion transversal del elemento    ');
        barras{cont}.altura = input('digite  la altura de la seccion transversal del elemento    ');
        barras{cont}.I=((barras{cont}.base*(barras{cont}.altura)^3)/12);
        % este condicional funciona para la formula de J recordar que en
        % esta formula t es el valor mas pequeño y b es el mayor esto
        % cambia el resultado de la formula J
        if barras{cont}.base >= barras{cont}.altura;   
        b = barras{cont}.base;
        t = barras{cont}.altura;
        barras{cont}.J=(1/3-0.21*((t/b)*(1-(1/12)*(t/b)^4)))*b*(t)^3;
        elseif barras{cont}.altura > barras{cont}.base;
        t=barras{cont}.base;
        b=barras{cont}.altura;
        barras{cont}.J=((1/3)-0.21*(t/b)*(1-(1/12)*(t/b)^4))*b*(t)^3;    
        end
        %disp (barras{cont}.J);% por si se necesita 
        %verificador si daba al J disp (barras{cont}.J) SI DA;
        disp (' ');      
     end
     fprintf ('PROPIEDADES DEL MATERIAL DEL ELEMENTO NUMERO %i.   ',cont );
     disp (' ');
     disp (' ');
     barras{cont}.E = input ('digite el modulo de elasticidad de este elemento (E)    '); 
     barras{cont}.G = input ('digite el modulo de cortante del elemento (G)    '); 
     disp (' ');     
    %matriz de rigidez del elemento
    % a continuacion estos if nos permiten saber que matriz de rigidez usar
    % recuerde que los elementos orientados en x-z se analisan con una matriz de rigidez diferente  
    % se mira si los elemntos son paralelos a los criterios de los planos
    % coordenados
    % criterio elemento orientado en x > vector a (1,0,0)
    % criterio elemento orientado en y > cector a (0,1,0)
    % se le solicitan las coordenadas de los dos nudos de la barra 
     disp (' ');
     fprintf ('COORDENADAS DEL NUDO INICIAL Y FINAL DE LA BARRA NUMERO %i.', cont );
     disp (' ');
     disp (' ');
     % en el ingreso de estas coordenadas recuerde que por facilidad se
     % deben ingresar los datos de las coordenadas  primero los mas
     % pequeños osea izquierda derecha y en y de abajo hacia arriba primero
     % las coordenadas mas pequeñas  a las mas grandes 
     % xin, yin son las coordenadas mas pequeñas
     disp ('recuerde que debe digitar las coordenadas primero el punto mas pequeño y despues el punto mas grande')
     disp (' ');
     barras{cont}.coordxin= input ('digite la coordenada en x inicial del elemento    ');
     barras{cont}.coordyin= input ('digite la coordenada en y inicial del elemento    ');
    % barras{cont}.coordzin= input ('digite la coordenada en z inicial del elemento    ');
     disp ('                 ');
     % xfin, yfin son las coordenadas mas grandes 
     barras{cont}.coordxfin= input ('digite la coordenada en x final del elemento    ');
     barras{cont}.coordyfin= input ('digite la coordenada en y final del elemento    '); 
    % barras{cont}.coordzfin= input ('digite la coordenada en z final del elemento    ');
     disp ('                 ');
     disp ('                 ');
    % vector posicion del elemento actual
   barras{cont}.vector =[barras{cont}.coordxfin-barras{cont}.coordxin,barras{cont}.coordyfin-barras{cont}.coordyin,0];
   barras{cont}.co = barras{cont}.vector(1)/norm(barras{cont}.vector); 
   barras{cont}.s = barras{cont}.vector(2)/norm(barras{cont}.vector);
   barras{cont}.MATRANSFORM = [barras{cont}.co -barras{cont}.s 0 0 0 0;
                            barras{cont}.s barras{cont}.co 0 0 0 0;
                            0 0 1 0 0 0;
                            0 0 0 barras{cont}.co -barras{cont}.s 0;
                            0 0 0 barras{cont}.s barras{cont}.co 0;
                            0 0 0 0 0 1];       
    %disp (barras{cont}.vector) no se si usarlo o si sirva;
    % esta condici
    %--------------------------------------------------------------------------------------------------------------corregir solo se necesita una sola matriz que es la de elementos orientados aleatoriamente.
    %condición paralelismo eje x 
    % elemento orientado en x
    if ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))== 1 | ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))== -1 ;
        disp ('elemento paralelo al eje x ');     
         k1=(barras{cont}.J * barras{cont}.G)/(norm(barras{cont}.vector));
         k2=(4*barras{cont}.E*barras{cont}.I)/(norm(barras{cont}.vector));
         k3=(6* barras{cont}.E * barras{cont}.I)/(norm(barras{cont}.vector))^2;
         k4=(12*barras{cont}.E * barras{cont}.I)/(norm(barras{cont}.vector))^3;
         % se quito el sub indice x barras{cont}.Krigidezx
         barras{cont}.Krigidez=[k1 0 0 -k1 0 0;
                                 0 k2 -k3 0 k5 k3;
                                 0 -k3 k4 0 -k3 -k4;
                                 -k1 0 0 k1 0 0;
                                 0 k5 -k3 0 k2 k3;
                                 0 k3 -k4 0 k3 k4];
    % condición paralelismio con el eje y 
    % elemento orientado en y
    elseif ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector))==1 | ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector))==-1 ; 
        disp ('elemento paralelo al eje y ');  
        k1=(barras{cont}.J*barras{cont}.G)/(norm(barras{cont}.vector));
        k2=(4*barras{cont}.E*barras{cont}.I)/(norm(barras{cont}.vector));
        k3=(6* barras{cont}.E * barras{cont}.I)/(norm(barras{cont}.vector))^2;
        k4=(12*barras{cont}.E * barras{cont}.I)/(norm(barras{cont}.vector))^3;
        k5=(2*barras{cont}.E * barras{cont}.I)/(norm(barras{cont}.vector));
        % se quito el sub indice y barras{cont}.Krigidezy 
        barras{cont}.Krigidez=[k2 0 k3 k5 0 -k3;
                                0 k1 0 0 -k1 0; 
                                k3 0 k4 k3 0 -k4;
                                k5 0 k3 k2 0 -k3;
                                0 -k1 0 0 k1 0;
                                -k3 0 -k4 -k3 0 k4];
    % si la barra no es paralela a ningun eje    
    elseif ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector)) ~= 1 | ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))~= 1 | ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector)) ~= -1 | ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))~= -1;
        barras{cont}.vector;
        disp ('este elemento no es paralelo ni a x ni a y');
        disp(' ');
        co = barras{cont}.vector(1)/norm(barras{cont}.vector); 
        s = barras{cont}.vector(2)/norm(barras{cont}.vector);
        %verificando si estaba tomando los valores correctos
        %disp ( barras{cont}.vector(1));
        %disp ( barras{cont}.vector(2));
        %disp (norm(barras{cont}.vector));
        %disp (co);
        %disp (s);
             MATRANSFORM = [co -s 0 0 0 0;
                            s co 0 0 0 0;
                            0 0 1 0 0 0;
                            0 0 0 co -s 0;
                            0 0 0 s co 0;
                            0 0 0 0 0 1];      
        % para facilitar la creacion de la matriz de  rigidez hago los
        % valores mas comunes como GJ/L Y EI/L COMO tor,flec.     
        tor = (barras{cont}.G * barras{cont}.J)/norm(barras{cont}.vector);     
        flec1 =(barras{cont}.E*barras{cont}.I)/norm(barras{cont}.vector); 
        flec2 = (barras{cont}.E*barras{cont}.I)/(norm(barras{cont}.vector))^2;
        flec3 = (barras{cont}.E*barras{cont}.I)/(norm(barras{cont}.vector))^3;
        if  barras{cont}.vector(2)/ barras{cont}.vector(2) == -1
           s = -s; 
           disp (s);
        end
        k1 = tor*co^2+4*flec1*s^2;
        k2 = co*s*(tor-4*flec1);
        k3 = 6*flec2*s;
        k4 = -tor*co^2 + 2*flec1*s^2;
        k5 = (tor+2*flec1)*co*s;       
        k6 = tor*s^2 + 4*flec1*co^2;
        k7 = 6*flec2*co;
        k8 = (tor+2*flec1)*co*s;
        k9 = -tor*s^2+2*flec1*co^2;
        k10 = 12*flec3;
        barras{cont}.Krigidez = [k1 k2 k3 k4 -k5 -k3;
                                 k2 k6 -k7 -k8 k9 k7;
                                 k3 -k7 k10 k3 -k7 -k10;
                                 k4 -k5 k3 k1 k2 -k3;
                                 -k8 k9 -k7 k2 k6 k7;
                                 -k3 k7 -k10 -k3 k7 k10];       
    end
    % en la siguiente parte del código se le pregunta al usuario que cargas
    % hay dentro de la viga y cuantas asumiendo que las cargas SIEMPRE son
    % en sentido negativo y recordando el nuevo cambio del sistema
    % coordenado los momento horarios son positivos y los antihorarios
    % negativos
    disp (' ');
    fprintf ('FUERZAS Y MOMENTOS QUE ACTUAN EN EL INTERIOR DE LA VIGA NUMERO %i.',cont);
    disp (' ');
    disp (' ');
    b = input('¿ Existen fuerzas distribuidas en este elemento ? [Y][N]:    ','s');
    c = input('¿ Existen fuerzas puntuales en este elemento ? [Y][N]:    ','s');
    w = input('¿ Existen momentos al interior de este elemento ? [Y][N]:    ','s');
    % a continuacion definimos la matriz de transformación del
    % elemento orientado aleatoriamente.
    %empiesa en cero este vector
    barras{cont}.accionesfijas = zeros(6,1);
    % se inican los valores de los momentos y las reacciones en 0 antes de
    % que se le pregunte al usuario el valor de estas cargas
    barras{cont}.momentocargadisizq = 0;
    barras{cont}.momentocargadisder = 0;
    barras{cont}.reaccioncargadisizq = 0;
    barras{cont}.reaccioncargadisder = 0;
    valorL = norm(barras{cont}.vector); 
    %para cargas distribuidas
    if b == 'y' | b == 'Y';  
        disp (' ');
        fprintf ('¿ Cuantas fuerzas distribuidas hay en el elemento %i ?.   ', cont);
        numfuerzdis = input ('  ');
        for d = 1 :  numfuerzdis; 
            disp (' ');
            fprintf ('TIPO DE CARGA DISTRIBUIDA NUMERO %i EN EL ELEMENTO NUMERO %i.', d , cont);
            disp (' ');
            disp (' ');
            disp ('1) carga rectangular ');
            disp (' ');
            e = input ('seleccione un opción    ');     
            % en el siguiente if lo que hago es dar la posibilidad de
            % ingresar las ecuaciones de una viga bajo distintas hipotesis
            % de carga.
            if e == 1 ;
              disp (' ');  
              fprintf ('CARGA DISTRIBUIDA NÚMERO %i APLICADA EN EL ELEMENTO NUMERO %i.', d , cont);
              disp (' ');
              disp (' ');
              disp (' ');
              valorcarga = input ('Digite el valor de la carga distribuida.    ' );
              distinicialdis = input ('Digite la distancia en la que empieza a aplicarse la carga distribuida.    ');
              distfinaldis = input ('Digite la distancia en la que termina de aplicarse la carga.    ');
              valorD  = distfinaldis - distinicialdis ;
              valorC = valorD /2;
              valordeA = valorC  +  distinicialdis;
              valordeB = valorL  - valordeA;        
              barras{cont}.momentocargadisizq = barras{cont}.momentocargadisizq-2*valorC*valorcarga*(3*valordeA*valordeB^2-valorC^2*(3*valordeB-1*valorL))/(3*(valorL)^2);   
              barras{cont}.momentocargadisder =  barras{cont}.momentocargadisder + 2*valorC*valorcarga*(3*valordeA^2*valordeB-3*valordeA*valorC^2+valorC^2*valorL)/(3*(valorL)^2);   
              barras{cont}.reaccioncargadisizq = barras{cont}.reaccioncargadisizq + 2*valorC*valorcarga*(2*valordeA^3-3*valordeA^2*valorL + 2*valordeA*valorC^2-(valorC^2-valorL^2)*valorL)/(valorL^3);           
              barras{cont}.reaccioncargadisder = barras{cont}.reaccioncargadisder -2*valorC*valorcarga*(2*valordeA^3-3*valordeA^2*valorL+2*valordeA*valorC^2-valorC^2*valorL)/(valorL^3); 
              % elemento paralelo a y.
              if ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector))==1 | ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector))==-1 ;
                  barras{cont}.accionesfijas = [-barras{cont}.momentocargadisizq;0;barras{cont}.reaccioncargadisizq;-barras{cont}.momentocargadisder;0;barras{cont}.reaccioncargadisder];
              % elemento paralelo a x.
            elseif ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))== 1 | ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))== -1 ;
                 barras{cont}.accionesfijas = [0; barras{cont}.momentocargadisizq;barras{cont}.reaccioncargadisizq;0; barras{cont}.momentocargadisder;barras{cont}.reaccioncargadisder];          
              %-----------------condicion barra sin paalelismo
              % fuerzas distribuidas
              elseif ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector)) ~= 1 | ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))~= 1 | ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector)) ~= -1 | ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))~= -1;                     
                 barras{cont}.accionesfijas = [0;barras{cont}.momentocargadisizq;barras{cont}.reaccioncargadisizq;0;barras{cont}.momentocargadisder;barras{cont}.reaccioncargadisder];
                 %disp ( barras{cont}.accionesfijas);
                 barras{cont}.accionesfijas = MATRANSFORM*barras{cont}.accionesfijas;                               
              end                            
            end                        
        end
        disp (' ');
        fprintf ('ACCIONES FIJAS DEL ELEMENTO NUMERO %i.', cont);              
        disp (' ');              
        disp (' ');
        disp (barras{cont}.accionesfijas);             
        disp (' ');
    end  
    % para cargas puntuales  
    if c == 'y' | c == 'Y';
        disp (' ');
        fprintf ('¿ Cuantas fuerzas puntuales hay en el elemento %i ?.   ', cont);
        numfuerzpunt = input ('  ');
        for f = 1:numfuerzpunt;
            disp (' ');
            fprintf ('CARGA PUNTUAL NÚMERO %i APLICADA EN EL ELEMENTO NÚMERO %i.', f , cont);
            disp (' ');
            disp (' ');
            disp (' ');
            valorcargapunt = input ('Digite el valor de la carga puntual.    ' );
            distizqpunt = input ('Digite la distancia de izquierda a derecha a la que se localiza la carga puntual.    ');
            valorApunt = distizqpunt;
            valorBpunt = valorL - valorApunt ;
            barras{cont}.momentocargadisizq = barras{cont}.momentocargadisizq - (valorcargapunt*valorApunt*valorBpunt^2)/(valorL^2);         
            barras{cont}.momentocargadisder = barras{cont}.momentocargadisder + (valorcargapunt*valorApunt^2*valorBpunt)/(valorL^2);           
            barras{cont}.reaccioncargadisizq = barras{cont}.reaccioncargadisizq + (-valorcargapunt*valorBpunt^2*(2*valorBpunt-3*valorL))/(valorL^3) ;           
            barras{cont}.reaccioncargadisder = barras{cont}.reaccioncargadisder + (-valorcargapunt*valorApunt^2*(2*valorApunt-3*valorL))/(valorL^3) ;
            % vector acciones fijas elemeneto paralelo a y
            if ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector))==1 | ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector))==-1 ;
                barras{cont}.accionesfijas = [-barras{cont}.momentocargadisizq;0;barras{cont}.reaccioncargadisizq;-barras{cont}.momentocargadisder;0;barras{cont}.reaccioncargadisder];                                                           
            % vector acciones fijas elemento paralelo a x. 
            elseif ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))== 1 | ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))== -1 ;                  
                barras{cont}.accionesfijas = [0; barras{cont}.momentocargadisizq;barras{cont}.reaccioncargadisizq;0; barras{cont}.momentocargadisder;barras{cont}.reaccioncargadisder];                             
            %-----------------condicion barra sin paralelismo es diagonal
            %fuerzas puntuales 
            % vector  acciones fijas            
            elseif ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector)) ~= 1 | ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))~= 1 | ((dot ([0,1,0],barras{cont}.vector))/norm(barras{cont}.vector)) ~= -1 | ((dot ([1,0,0],barras{cont}.vector))/norm(barras{cont}.vector))~= -1;                 
               disp ('');
               barras{cont}.accionesfijas = [0;barras{cont}.momentocargadisizq;barras{cont}.reaccioncargadisizq;0;barras{cont}.momentocargadisder;barras{cont}.reaccioncargadisder];
               %disp ( barras{cont}.accionesfijas);
               barras{cont}.accionesfijas = MATRANSFORM*barras{cont}.accionesfijas;                                       
            end                                            
        end        
              disp (' ');
              fprintf ('ACCIONES FIJAS DEL ELEMENTO NÚMERO %i.', cont);              
              disp (' ');              
              disp (' ');
              disp (barras{cont}.accionesfijas);             
              disp (' ');                
    end    
    % en caso de que hallan momentos dentro de la viga    
    if w == 'Y' | w == 'y';        
        disp (' ');
        disp ('Por favor digite sobre cual eje del elemento se esta aplicando este momento');
        o = input ('indique sobre que eje por favor (x) o (y).    ','s');        
        % en dado caso que sea sobre el eje x del elemento  este eje es
        % local entonces sera un momento torsor---------------------------------------------------------NO ESTA COMPLETO PARA UN MOMENTO TORSOR SOBRE EL EJE DEL ELEMENTO                       
        if o == 'x' | o == 'X';        
        end        
        % si el momento es sobre el eje y local del elemento este sera un
        % momento normal
        % ---------------------------------------------------------------------------------------                      
        if o == 'Y' | o == 'y';           
        disp (' ');
        fprintf ('¿ Cuantas momentos sobre el eje Y local hay en el elemento %i ?.   ', cont);
        nummoment = input ('  ');            
        for r = 1 : nummoment;                        
            disp (' ');
            fprintf ('MOMENTO NUMERO %i APLICADO EN EL ELEMENTO NUMERO %i.', r , cont);
            disp (' ');
            disp (' ');            
            disp (' ');
            valormoment = input ('Digite el valor del momento.    ' );
            dismomenta = input ('Digite la distancia de izquierda a derecha donde se aplica el momento.    ');
            dismomentb = norm(barras{cont}.vector) - dismomenta;            
            barras{cont}.momentocargadisizq = barras{cont}.momentocargadisizq + valormoment*(dismomentb/norm(barras{cont}.vector))*(2-3*(dismomentb/norm(barras{cont}.vector)));                        
            barras{cont}.momentocargadisder = barras{cont}.momentocargadisder + valormoment*(dismomenta/norm(barras{cont}.vector))*(2-3*(dismomenta/norm(barras{cont}.vector)));                        
            barras{cont}.reaccioncargadisizq = barras{cont}.reaccioncargadisizq + valormoment*((6*dismomenta*dismomentb)/(norm(barras{cont}.vector))^3);                        
            barras{cont}.reaccioncargadisder = barras{cont}.reaccioncargadisder -valormoment*((6*dismomenta*dismomentb)/(norm(barras{cont}.vector))^3);            
            % sale mucho mejor trabajar todas las acciones fijas bajo un
            % sistema local de corrdenadas y SIEMPRE multiplicarle la
            % matriz de transformacion de un sistema local a global que en
            % este caso es barras{cont}.MATRANSFORM.            
            barras{cont}.accionesfijas = [0; -barras{cont}.momentocargadisizq; barras{cont}.reaccioncargadisizq; 0; -barras{cont}.momentocargadisder; barras{cont}.reaccioncargadisder];
            barras{cont}.accionesfijas = barras{cont}.MATRANSFORM * barras{cont}.accionesfijas;                       
        end        
              disp (' ');
              fprintf ('ACCIONES FIJAS DEL ELEMENTO NUMERO %i.', cont);              
              disp (' ');              
              disp (' ');
              disp (barras{cont}.accionesfijas);             
              disp (' ');                  
        end          
    end        
    % el siguiente for le permite al usuario la ubicación de los grados de
    % libertad de la estructura de izquierda a derecha de menor a mayor    
    disp ('');
    fprintf ('ASIGNACION DE LOS GRADOS DE LIBERTAD DEL ELEMENTO NUMERO %i:',cont);
    disp (' ');
    disp (' ');
    disp ('por favor digite el valor del grado de libertad de cada nudo del elemento en el siguiente orden:');
    disp ('Mxi, Myi, Zi, Mxj, Myj, Zj');
    disp (' ');    
    for x = 1 : 6;        
        barras{cont}.gradosdelibertad(x)=input('');        
    end 
end 
% se ensambla la matriz de acciones fijas de la estructura. 
F=zeros(nodos*3,1); 
for cont = 1 : barrasnum;    
    barras{cont}.NUMERODDEGDL = size(barras{cont}.accionesfijas,1);    
    for i=1:barras{cont}.NUMERODDEGDL        
        ifila=barras{cont}.gradosdelibertad(i);
        F(ifila,1)=F(ifila,1)+ barras{cont}.accionesfijas(i);    
    end   
end 
% se pregunta las fuerzas que actuan en los nodos de la estructura
disp (' ');
numf = input ('¿ Cuantas fuerzas actuan en los nodos de la estructura ?.    ');
fuerznodos = zeros (nodos*3-dc,1); 
for x = 1 :numf    
    disp (' ');
    fprintf ('Indique el grado de libertad donde se ubica la fuerza numero %i.',x);
    disp (' ');
    disp (' ');
    psfuervec = input ('');    
    disp (' ');
    fprintf ('Ingrese el valor de la fuerza que actúa en el grado de libertad numero %i.    ',psfuervec); 
    disp (' ');
    disp (' ');
    funu = input (' ');    
    fuerznodos(psfuervec)= funu;    
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
         DB(posivec-(nodos*3-dc))=desplaz;
    end
else
    DB=zeros(dc,1);
end 
% En el siguiente ciclo for se ensambla la matriz de rigidez de la estructura 
KGENERAL=zeros(nodos*3,nodos*3); 
for cont = 1 : barrasnum;   
    barras{cont}.NUMERODDEGDL = size(barras{cont}.Krigidez);   
    for i=1 : barras{cont}.NUMERODDEGDL;       
        ifil = barras{cont}.gradosdelibertad(i);        
        for j=1 : barras{cont}.NUMERODDEGDL;        
            icol = barras{cont}.gradosdelibertad(j);            
            KGENERAL(ifil , icol) = KGENERAL(ifil, icol) + barras{cont}.Krigidez(i,j) ;            
        end        
    end    
end 
% se copia la parte RAA de la matriz de rigidezes. 
RAA =zeros (nodos*3-dc,nodos*3-dc); 
for i = 1 : nodos*3-dc    
    for j =1 : nodos*3-dc  
        RAA (i,j) = KGENERAL (i,j);
    end   
end 
%se copia la parte RAB  de la matriz de rigidezes. 
RAB =zeros(nodos*3-dc,dc); 
x = 1; % contadores 
y = 1;
for i = 1:nodos*3-dc
    for j = (nodos*3-dc + 1):nodos*3
        RAB(x,y) = KGENERAL (i,j);  
    end
    x=x+1;
    y=1;
end
% se copia la parte de RBA de la matriz de rigidezes.
RBA = zeros (dc,nodos*3-dc);
x = 1;
y = 1;
for i = (nodos*3-dc + 1):nodos*3
    for j = 1 : (nodos*3-dc)
        RBA (x,y) = KGENERAL(i,j); 
        y = y+1;
    end
    y = 1;
    x= x+1; 
end
% se copia la parte de RBB de la matriz de rigidezes.
RBB = zeros(dc,dc);
x = 1;
y = 1;
    for i=(nodos*3-dc+1):nodos*3
        for j = (nodos*3-dc+1):nodos*3
        RBB (x,y) = KGENERAL(i,j);     
        y = y + 1;
        end
        x = x + 1;
        y =1; 
    end 
% seguidamente se hace lo mismo con el vector de acciones fijas . 
% vector de acciones fijas Aa0.
Aa0=zeros(nodos*3-dc,1);
     for i=1:nodos*3-dc
        Aa0(i)=F(i);
    end 
% vector de acciones fijas Ab0 
Ab0 = zeros (dc,1);
x = 1;
    for i = (nodos*3-dc+1):nodos*3
        Ab0(x) = F(i); 
        x = x + 1;
    end
% se plantean y se resuelven las ecuaciones para encontrar los
% desplazamientos desconocidos y las fuerzas desconocidas
deltasdesconocidos = (inv(RAA))*(fuerznodos-Aa0-RAB*DB);
ABN = Ab0 + RBA*deltasdesconocidos  + RBB*DB; 
% se construye el vector desplazamienos completo------------------nuevo
% linea
y=1;
z=1;
Ddefcom=zeros(nodos*3,1);
    for x=1:(nodos*3)
        if(x<=(nodos*3-dc))
            Ddefcom(x)=deltasdesconocidos(y);
            y=y+1;
        else
            Ddefcom(x)=DB(z);
            z=z+1;
        end
    end
% Se hallan las acciones de las fuerzas internas en cada elemento de la
% estructura.
    for cont = 1 : barrasnum;
        barras{cont}.desplas=zeros(6,1);
    end
    for cont = 1 : barrasnum;
        for x = 1 : 6;
            i=barras{cont}.gradosdelibertad(x);
            barras{cont}.desplas(x)=Ddefcom(i); 
        end
        barras{cont}.Fintgen = barras{cont}.Krigidez*barras{cont}.desplas+barras{cont}.accionesfijas;
    % no sirve esto ----- barras{cont}.Fintloc = barras{cont}.Krigidez*barras{cont}.desplas+barras{cont}.accionfijalocal;
    end
% matriz de accines en cada elemento localmente habalndo
        for cont = 1 : barrasnum; 
            barras{cont}.Fintloc=zeros(6,1);
        end   
        for cont = 1: barrasnum;
            barras{cont}.Fintloc=(inv(barras{cont}.MATRANSFORM))* barras{cont}.Fintgen; 
        end
% en el siguiente switch menu se crea un menu para el usuario
h = 1;
while h == 1;
g = menu ('OPCIONES','salir.','Matriz de rigidez de cada elemento.','Matriz de rigides de la estructura.','Vector de acciones fijas de la estructura.','Desplazamientos desconocidos.','Fuersas desconocidas.','vector acciones en los nodos de la estructura','fuerzas internas en cada elemento coordenadas generlaes','Fuerzas internas en cada elemento coordenadas locales' );
        switch g ;
        case 1
            break
        case 2
            disp (' ');
            disp ('MATRIZ DE RIGIDEZ DE CADA ELEMENTO:')
            disp (' ');
            disp ('Elementos orientados en x y en y');
%en el siguiente for se muestra la matriz de rigidez de cada elemento 
        for cont = 1 :1:barrasnum
            disp (' ');
            fprintf ('Elemento numero %i ', cont);
            disp (' ');
            disp (' ');
            disp ((barras{cont}.Krigidez));
            disp (' ');
        end
        case 3
            disp (' ');
            disp ('MATRIZ DE RIGIDEZ DE LA ESTRUCTURA ENSAMBLADA COMPLETAMENTE:');
            disp (' ');
            % matriz de rigidez de la estructura
            disp (' ');
            disp (KGENERAL); 
        case 4    
            disp (' ');
            disp ('VECTRO DE ACCIONES FIJAS DE LA ESTRUCTURA');
            disp (' ');
            disp (' ');
            disp (F);
        case 5
            disp (' ');
            disp ('VECTOR DESPLAZAMIENTOS DESCONOCIDOS DE LA ESTRUCTURA');
            disp (' ');
            disp (' '); 
            disp( deltasdesconocidos );
        case 6
            disp (' ');
            disp ('VECTOR FUERSAS DESCONOCIDAS DE LA ESTRUCTURA');
            disp (' ');
            disp (' '); 
            disp( ABN );
        case 7
            disp (' ');
            disp ('VECTOR ACCIONES EN LOS NODOS DE LA ESTRUCTURA');
            disp (' ');
            disp (' ');     
            disp(fuerznodos);    
        case 8
            disp (' '); 
            disp ('FUERZAS INTERNAS EN CADA ELEMENTO COORDENADAS GENERALES');
            disp (' ');
            for cont = 1: barrasnum 
                disp (' ');
                fprintf('fuerzas internas elemento numero %i.',cont);
                disp (' ');
                disp (' ');
                disp(barras{cont}.Fintgen);  
            end
        case 9
            disp (' '); 
            disp ('FUERZAS INTERNAS EN CADA ELEMENTO COORDENADAS LOCALES');
            disp (' ');
            for cont=1: barrasnum;
                disp (' ');
                fprintf('fuerzas internas elemento numero %i.',cont);
                disp (' ');
                disp (' ');
                disp(barras{cont}.Fintloc);     
            end          
        end
    end
end

