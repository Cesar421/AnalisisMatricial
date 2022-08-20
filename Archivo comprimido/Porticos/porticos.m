

%CALCULO DE REACCINES PORTICOS
%este programa analiza PORTICOS con acciones exclusivamente sobre el eje y respecto a un plano
%coordenado global.


%clc;
 
disp ('Recuerde manejar siempre las mismas unidades en el ingreso de datos');
disp ('Ademas recuerde ingresar las coordenadas de los nodos SIEMPRE de izquierda a derecha del menor a mayor punto');
disp ('cuando ingrese los grados de libertad se empiesa siempre a enumerar de acuerdo a los desplazamientos desconocidos')
disp ('  ');
disp ('  ');
 
disp ('                                 INICIO        ');
disp ('Programa que calcula reacciones en un PORTICO ');
c=input ('digite el número de rotulas   '); 
barrasnum =input ('digite el número de barras    ');
nodos =input ('digite el número de nodos    ');
dc =input ('digite el número de desplazamientos conocidos    ');
dc=dc+c
GIE=barrasnum+dc-3*nodos-c;
GIC = 3*nodos-dc;
disp (' ');
disp ('grado de indeterminación cinematico    ');
fprintf ('GIC = %i.  ',GIC); 
disp ('grado de indeterminación estatico   ');
fprintf ('GIC = %i.  ',GIE); 
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
        barras{cont}.L=input('digite la longitud del elemento:    ')
        disp (' ');
        
     %seccion rectangular
     
     elseif a == 2
         
        disp (' ');
        IO=input('si conoce las proiedades digite 1:    ')
       
        if IO==1
        barras{cont}.A=input('digite el Area    ')
        if barras{cont}.A==0
        barras{cont}.A = input ('digite el Area:   '); 
        end  
        barras{cont}.I=input('digite la inercia:  ');
        if barras{cont}.I==0
        barras{cont}.I = input ('digite lainercia:  '); 
        end
        barras{cont}.E =input ('digite el modulo de elasticidad de este elemento (E) :   '); 
        if barras{cont}.E==0
        barras{cont}.E = input ('digite el modulo de elasticidad de este elemento (E) :    '); 
        end 
        else
        barras{cont}.base=.4%input('digite la BASE de la seccion transversal del elemento:    ');
        barras{cont}.altura =.6%input('digite  la AlTURA(H) de la seccion transversal del elemento:  ');
        barras{cont}.I=((barras{cont}.base*(barras{cont}.altura)^3)/12);
        barras{cont}.A = barras{cont}.base*barras{cont}.altura;
        barras{cont}.E =input ('digite el modulo de elasticidad de este elemento (E):   '); 
        if barras{cont}.E==0
        barras{cont}.E = input ('digite el modulo de elasticidad de este elemento (E):   '); 
        
        end
        
        if barras{cont}.base >= barras{cont}.altura;
            
            
        b = barras{cont}.base;
        t = barras{cont}.altura;
        barras{cont}.J=(1/3-0.21*((t/b)*(1-(1/12)*(t/b)^4)))*b*(t)^3;
        barras{cont}.area=b*t;
        
        elseif barras{cont}.altura > barras{cont}.base;
            
        t=barras{cont}.base;
        b=barras{cont}.altura;
        barras{cont}.J=((1/3)-0.21*(t/b)*(1-(1/12)*(t/b)^4))*b*(t)^3;
       
        end
        end
        barras{cont}.L=input('digite la longitud del elemento:   ')
        % este condicional funciona para la formula de J recordar que en
        % esta formula t es el valor mas pequeño y b es el mayor esto
        % cambia el resultado de la formula J
        %disp (barras{cont}.J);% por si se necesita 
        %verificador si daba al J disp (barras{cont}.J) SI DA;
        
        disp (' ');
                
     end
     
     
     fprintf ('PROPIEDADES DEL MATERIAL DEL ELEMENTO NUMERO %i.   ',cont );
     disp (' ');
     disp (' ');
     fprintf ('BARRA %i E = %i.    ',cont,barras{cont}.E); 
     fprintf ('BARRA %i I = %i.    ',cont,barras{cont}.I);
     disp ('  ');
     disp ('  ');
     
     fprintf ('BARRA %i L = %i.    ',cont,barras{cont}.L);
     fprintf ('BARRA %i A = %i.    ',cont,barras{cont}.A);

     disp (' '); 
     angulo=input ('digite el ANGULO por eje local al global en sentido antihorario :    ')
     rad=(angulo*pi()/180)
     barras{cont}.cos=cos(rad);
     barras{cont}.sen=sin(rad);
     C=barras{cont}.cos
     S=barras{cont}.sen
     Axial=(barras{cont}.A*barras{cont}.E)/(barras{cont}.L);
     Flector14=(4*(barras{cont}.I*barras{cont}.E))/(barras{cont}.L);
     Flector2=(6*(barras{cont}.I*barras{cont}.E))/(barras{cont}.L)^2;
     Flector3=(12*(barras{cont}.I*barras{cont}.E))/(barras{cont}.L)^3;
     Flector12=(2*(barras{cont}.I*barras{cont}.E))/(barras{cont}.L);
     Flector33=(3*(barras{cont}.I*barras{cont}.E))/(barras{cont}.L)^3;
     Flector32=(3*(barras{cont}.I*barras{cont}.E))/(barras{cont}.L)^2;
     Flector31=(3*(barras{cont}.I*barras{cont}.E))/(barras{cont}.L);
     disp('es este un elemento con ROTULA?')
     disp ('Al INICIO Digite 1')
     disp ('Al FINAL Digite 2')
     disp ('Al en AMBOS extremos  Digite 3')
     disp ('ELEMENTO RIGIDO Y/O FLEXIBELE 4')
     condicionalrotula=input('digite su respueta   ')
     if condicionalrotula==1 
     barras{cont}.Krigidez=[Axial 0 0 -Axial 0 0;
                            0 Flector33 0 0 -Flector33 Flector32;
                            0 0 0 0 0 0;
                            -Axial 0 0 Axial  0 0;
                            0 -Flector33 0 0 Flector33 -Flector32;
                            0 Flector32 0 0 -Flector32 Flector31]
     
     elseif condicionalrotula==2
              barras{cont}.Krigidez=[Axial 0 0 -Axial 0 0;
                            0 Flector33 Flector32  0 -Flector33 0;
                            0 Flector32 Flector31 0 -Flector32 0;
                            -Axial 0 0 Axial  0 0;
                            0 -Flector33 -Flector32 0 Flector33 0;
                            0 0 0 0 0 0]
     
         elseif condicionalrotula==3
          barras{cont}.Krigidez=[Axial 0 0 -Axial 0 0;
                            0 0 0 0 0 0;
                            0 0 0 0 0 0;
                            -Axial 0 0 Axial  0 0;
                            0 0 0 0 0 0;
                            0 0 0 0 0 0]
                        
         
         
    else        
     barras{cont}.Krigidez=[(Axial*C^2+Flector3*S^2) (Axial-Flector3)*C*S      -Flector2*S  -(Axial*C^2+Flector3*S^2) -(Axial-Flector3)*C*S  -Flector2*S;
                           (Axial-Flector3)*C*S   Axial*S^2+Flector3*C^2      Flector2*C  -(Axial-Flector3)*C*S  -(Axial*S^2+Flector3*C^2) Flector2*C;
                           -Flector2*S                -Flector2*C              Flector14      Flector2*S              -Flector2*C             Flector12;
                           -(Axial*C^2+Flector3*S^2) -(Axial-Flector3)*C*S      Flector2*S    Axial*C^2+Flector3*S^2   (Axial-Flector3)*C*S  Flector2*S;
                           -(Axial-Flector3)*C*S    -(Axial*S^2+Flector3*C^2) -Flector2*C     (Axial-Flector3)*C*S     Axial*S^2+Flector3*C^2  -Flector2*C;                
                           -Flector2*S                  Flector2*C            Flector12       Flector2*S                 -Flector2*C           Flector14 ]                                             
     K=barras{cont}.Krigidez
     end 
     if condicionalrotula==4
         alphai=input('digite el valor ALPHA i'  )
         alphaj=input('digite el valor ALPHA j'  )
         alpha=alphai*alphaj+4*alphai+4*alphaj+12
         d1=(alphai*alphaj+alphai+alphaj)/alpha
         d2=(alphai*alphaj+2*alphai)/alpha
         d3=(alphai*alphaj+2*alphaj)/alpha
         d4=alphai*alphaj/alpha
         d5=(alphai*alphaj+3*alphai)/alpha 
         d6=(alphai*alphaj+3*alphaj)/alpha
         barras{cont}.Krigidez=[Axial  0             0    -Axial    0          0;
                                0   d1*Flector3 d2*Flector2 0 -d1*Flector3 d3*Flector2;
                                0   d2*Flector2 d5*Flector14 0 -d2*Flector2 d4*Flector12;
                                -Axial  0     0    Axial    0          0;
                                0 -d1*Flector3 -d2*Flector2 0 d1*Flector3 -d3*Flector2;
                                0 d3*Flector2 d4*Flector12 0 -d3*Flector2 d6*Flector14;]
                            %%
                            % 
                            % $$e^{\pi i} + 1 = 0$$
                            % 
                            % $$e^{\pi i} + 1 = 0$$
                            % 
                            %  PREFORMATTED
                            %  TEXT
                            % 
                            % 
                            % 
     
     
     end
  
     barras{cont}.T=[C -S 0 0 0 0;
                     S C 0 0 0 0;
                     0 0 1 0 0 0;
                     0 0 0 C -S 0;
                     0 0 0 S C 0;
                     0 0 0 0 0 1];
     
     TRANSFORMACION=barras{cont}.T
     if condicionalrotula==1 | condicionalrotula==2 | condicionalrotula==3 | condicionalrotula==4
     barras{cont}.Krigidez=TRANSFORMACION*barras{cont}.Krigidez*inv(TRANSFORMACION);
     K=barras{cont}.Krigidez
     end
    %%Fuerzas fijas delelemento i
    
    disp (' ');
    fprintf ('FUERZAS Y MOMENTOS QUE ACTUAN EN EL INTERIOR DE LA VIGA NUMERO %i.',cont);
    disp (' ');
    disp (' ');
    b = input('¿ Existen fuerzas distribuidas en este elemento ? [Y][N]:    ','s');
    c = input('¿ Existen fuerzas puntuales en este elemento ? [Y][N]:    ','s');
    %empiesa en cero este vector
    
    barras{cont}.accionesfijas = zeros(6,1);
    
    % se inican los valores de los momentos y las reacciones en 0 antes de
    % que se le pregunte al usuario el valor de estas cargas
    
    barras{cont}.momentoizq = 0;
    barras{cont}.momentoder = 0;
    barras{cont}.reaccionizq = 0;
    barras{cont}.reaccionder = 0;
    barras{cont}.xizq=0;
    barras{cont}.xder=0;
    %carga distribuida
    if b == 'y' | b == 'Y';  
        
        disp (' ');
        fprintf ('¿ Cuantas fuerzas distribuidas hay en el elemento %i ?.   ', cont);
        numfuerzdis = input ('  ');
       
        for d = 1 :  numfuerzdis; 
            
            disp (' ');
            fprintf ('TIPO DE CARGA DISTRIBUIDA NUMERO %i EN EL ELEMENTO NUMERO %i.', cont);
            disp (' ');
            disp (' ');
            disp ('CARGA DISTRIBUIDA RECTANGULAR ');
            disp (' ');
                   
              disp (' ');  
              fprintf ('CARGA DISTRIBUIDA NÚMERO %i APLICADA EN EL ELEMENTO NUMERO %i.', d , cont);
              disp (' ');
              disp (' ');
              
              disp (' ');
              w= input ('Digite el valor de la carga distribuida:  ' );
              disp ('  ');
    
              valorD  =input('Digite la distancia de aplicacion de la carga:  ');
              valorT=w*valorD;
              barras{cont}.momentoizq =(-w*valorD^2)/12;
              barras{cont}.momentoder = (w*valorD^2)/12;
              barras{cont}.reaccionizq = -valorT/2;
              barras{cont}.reaccionder = -valorT/2;            
            end
        end
   
    %carga puntual
    if c == 'y' | c == 'Y';           
         disp (' ');
        fprintf ('¿ Cuantas fuerzas puntuales hay en el elemento %i ?.   ', cont);
        nump = input ('  ');
        for d = 1 :  nump;
              fprintf ('carga puntual del elemento %i ?.   ', cont);
              p= input ('Digite el valor de la carga PUNTUAL:  ' );
              disp('los momentos se calcularan (P(a)(b^2))/L^2 donde a es la distancia del apoyo izquierdo y b del apoyo derecho ')
              a  =input('Digite la distancia "a" de aplicacion de la carga:  ');
              b =input('Digite la distancia "b" de aplicacion de la carga:  ');
              l=barras{cont}.L
               barras{cont}.momentoizq=barras{cont}.momentoizq-((p*a*b^2)/(l^2));
               barras{cont}.momentoder =barras{cont}.momentoder+ ((p*a*b^2)/(l^2));
               barras{cont}.reaccionizq =barras{cont}.reaccionizq-p/2;
               barras{cont}.reaccionder =barras{cont}.reaccionder-p/2;  
                       
        end
    end 
    if condicionalrotula==1
    barras{cont}.accionesfijas(1,1)=barras{cont}.xizq;     
    barras{cont}.accionesfijas(2,1)=barras{cont}.reaccionizq-((1.5*barras{cont}.momentoizq)/barras{cont}.L);
    barras{cont}.accionesfijas(3,1)=0;
    barras{cont}.accionesfijas(4,1)=barras{cont}.xder;
    barras{cont}.accionesfijas(5,1)=barras{cont}.reaccionder+((1.5*barras{cont}.momentoizq)/barras{cont}.L);
    barras{cont}.accionesfijas(6,1)=barras{cont}.momentoder-(0.5*barras{cont}.momentoizq);
    elseif condicionalrotula==2
    barras{cont}.accionesfijas(1,1)=barras{cont}.xizq;     
    barras{cont}.accionesfijas(2,1)=barras{cont}.reaccionizq-((1.5*barras{cont}.momentoder)/barras{cont}.L);
    barras{cont}.accionesfijas(3,1)=barras{cont}.momentoizq-(0.5*barras{cont}.momentoder);
    barras{cont}.accionesfijas(4,1)=barras{cont}.xder;
    barras{cont}.accionesfijas(5,1)=barras{cont}.reaccionder+((1.5*barras{cont}.momentoder)/barras{cont}.L);
    barras{cont}.accionesfijas(6,1)=0;
        
        elseif condicionalrotula==3
              barras{cont}.accionesfijas(1,1)=barras{cont}.xizq;     
    barras{cont}.accionesfijas(2,1)=barras{cont}.reaccionizq-((barras{cont}.momentoder+barras{cont}.momentoizq)/barras{cont}.L);
    barras{cont}.accionesfijas(3,1)=0;
    barras{cont}.accionesfijas(4,1)=barras{cont}.xder;
    barras{cont}.accionesfijas(5,1)=barras{cont}.reaccionder+(1.5*barras{cont}.momentoder/barras{cont}.L);
    barras{cont}.accionesfijas(6,1)=0;
        elseif condicionalrotula==4
    barras{cont}.accionesfijas(1,1)=barras{cont}.xizq;     
    barras{cont}.accionesfijas(2,1)=barras{cont}.reaccionizq;
    barras{cont}.accionesfijas(3,1)=barras{cont}.momentoizq;
    barras{cont}.accionesfijas(4,1)=barras{cont}.xder;
    barras{cont}.accionesfijas(5,1)=barras{cont}.reaccionder;
    barras{cont}.accionesfijas(6,1)=barras{cont}.momentoder; 
        
        else 
    barras{cont}.accionesfijas(1,1)=barras{cont}.xizq;     
    barras{cont}.accionesfijas(2,1)=barras{cont}.reaccionizq-(6/alpha*barras{cont}.L*(barras{cont}.momentoizq*(2+alphaj)+barras{cont}.momentoder*(2+alphai)));
    barras{cont}.accionesfijas(3,1)=(alphai/alpha)*(barras{cont}.momentoizq*(4+alphaj)-2*barras{cont}.momentoder);
    barras{cont}.accionesfijas(4,1)=barras{cont}.xder;
    barras{cont}.accionesfijas(5,1)=barras{cont}.reaccionder+(6/alpha*barras{cont}.L*(barras{cont}.momentoizq*(2+alphaj)+barras{cont}.momentoder*(2+alphai)));
    barras{cont}.accionesfijas(6,1)=(alphaj/alpha)*(barras{cont}.momentoder*(4+alphai)-2*barras{cont}.momentoizq);
    end  
    disp (' ')
      inclinado=input('es este un elemento inclinado digite 1 si asi lo es?    ')
      disp ('  ');
     disp ('  ');
       
      if inclinado==1
          barras{cont}.accionesfijas=barras{cont}.T*barras{cont}.accionesfijas
      
      else
          fprintf ('ACCIONES FIJAS DEL ELEMENTO NÚMERO %i.', cont)             
              disp (' ')
              disp (' ')
              barras{cont}.accionesfijas          
              disp (' ') 
              
      end 
      barras{cont}.accionesfijas
    %%--------------------------------------------------------------------------------
    disp ('');
    fprintf ('ASIGNACION DE LOS GRADOS DE LIBERTAD DEL ELEMENTO NUMERO %i:',cont);
    disp (' ');
    disp (' ');
    disp ('por favor digite el valor del grado de libertad de cada nudo del elemento en el siguiente orden:');
    disp ('Xi,Yi,Ti,Xj,Yj,Tj ');
    disp (' ');
    
    for x = 1 : 6;
        
        
        barras{cont}.gradosdelibertad(x)=input('');
    
    end
end
% se ensambla la matriz de acciones fijas de la estructura.
 

F=zeros(nodos*3,1);
 
for cont = 1 : barrasnum;
    cont
    barras{cont}.NUMERODDEGDL = size(barras{cont}.accionesfijas,1)
    
    for i=1:barras{cont}.NUMERODDEGDL
        
        ifila=barras{cont}.gradosdelibertad(i);
        F(ifila,1)=F(ifila,1)+ barras{cont}.accionesfijas(i);
     
    end
    F
end
%%fuerzas externas
    Fext=zeros(nodos*3-dc,1)
     disp (' ');
    numf = input ('¿ Cuantas fuerzas actuan en los nodos de la estructura ?.    ');
    Fext=zeros(nodos*3-dc,1)
 
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
Fext
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

KGENERAL=zeros(nodos*3,nodos*3);
 
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
        
    RBB (x,y) = KGENERAL(i,j)  ;   
        y = y + 1;
        
    end
    
    x = x + 1;
    y =1;
    
end
% seguidamente se hace lo mismo con el vector de acciones fijas . 
% vector de acciones fijas Aa0.
 
FA=zeros(nodos*3-dc,1);
 
for i=1:nodos*3-dc
    
     FA(i)=F(i);
 
end
 
% vector de acciones fijas Ab0
 
FB = zeros (dc,1);
x = 1;
 
for i = (nodos*3-dc+1):nodos*3
    
    FB(x) = F(i); 
    x = x + 1;
    
end
%%Fuerzas externas
FextA=zeros(nodos*3-dc,1);
 
for i=1:nodos*3-dc
    
     FextA(i)=Fext(i);
 
end
 
% vector de acciones fijas Ab0
 
%FextB = zeros (nodos*3-dc,1);

%for i = (nodos*3-dc+1):dc
    
 %   FextB(x) = Fext(i); 
  % x = x + 1;
    
%end

FextA
%FextB
FB
FA
RAA
RAB
RBA
RBB

% se plantean y se resuelven las ecuaciones para encontrar los
% desplazamientos desconocidos y las fuerzas desconocidas

deltasdesconocidos = (inv(RAA))*(FextA-FA)
 
ABN =  RBA*deltasdesconocidos  + (FB)

%%jj
y=1;
z=1;
Dcom=zeros(nodos*3,1);

for x=1:(nodos*3)
    
    if(x<=(nodos*3-dc))
        Dcom(x)=deltasdesconocidos(y);
        y=y+1;
        
    else
        
        Dcom(x)=DB(z);
        z=z+1;
        
    end
    
end
disp ('  ');
Dcom
disp ('  ');
     
% Se hallan las acciones de las fuerzas internas en cada elemento de la
% estructura.
for cont = 1 : barrasnum;
    
    barras{cont}.desplas=zeros(6,1);
    
end
disp('Fuerzas Globales ')
for cont = 1 : barrasnum;
    
    for x = 1 : 6;
        
       i=barras{cont}.gradosdelibertad(x);
       barras{cont}.desplas(x)=Dcom(i); 
        
    end
    
    barras{cont}.FintG = barras{cont}.Krigidez*barras{cont}.desplas+barras{cont}.accionesfijas;
    % no sirve esto ----- barras{cont}.Fintloc = barras{cont}.Krigidez*barras{cont}.desplas+barras{cont}.accionfijalocal;
    disp(cont)
    disp(barras{cont}.FintG)
end
disp('Fuerzas locales ')

% matriz de accines en cada elemento local

 for cont = 1 : barrasnum; 
     
     barras{cont}.Fintloc=zeros(6,1);
 end
    
for cont = 1: barrasnum;
    
    barras{cont}.FintL=(inv(barras{cont}.T))* barras{cont}.FintG;
    disp(cont)
    barras{cont}.FintL
end

%%-------------------------------------------------------------------
