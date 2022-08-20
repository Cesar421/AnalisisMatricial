%clc;
%clear;
%clear all;
%Datos iniciales de entrada
nnd=4;
nel=3;
nsec=3;
nmat=1;
GIC=6;
GL=3*nnd;
nrest=GL-GIC;
nFext=0;

%Dimensionamiento de matrices
Coord=zeros(nnd,2);
CAp=zeros(nnd,3);
Secc=zeros(nsec,2);
Mat=zeros(nmat,1);
Fext=zeros(nnd,3);
Elem=zeros(nel,5);
MGL=zeros(nnd,3);
KT=zeros(GL,GL);
GLel=zeros(nel,6);
Femporamiento=zeros(GL,1);

%Entrada de datos de los nudos
for i=1:1:nnd;
    nudo=i;
    Coord(i,1)=input(['Coordenada x del nudo ' num2str(i) '\n']);
    Coord(i,2)=input(['Coordenada y del nudo ' num2str(i) '\n']);  
end

%Entrada de datos de las restricciones
for ires=1:1:nrest;
    i=input(['Nudo restringido:' '\n']);
    j=input(['restringido en x=1, en y=2 o en giro=3?? ']);
    CAp(i,j)=1;
    disp (CAp);
end

%Matriz de grados de libertad
% 1 = apoyo restringido, 0 = apoyo sin restricción.
cont1=1;
cont2=GIC+1;
for i=1:1:nnd
    for j=1:1:3
        if CAp(i,j)== 1;
           MGL(i,j)=cont2;
           cont2=cont2+1;
        else
           MGL(i,j)=cont1;
           cont1=cont1+1;
        end
    end
end
disp (MGL);

%Fuerzas externas
for cont1=1:1:nFext
    i=input(['Nudo con carga' '\n']);
    j=input(['Dirección de la carga en x=1, y=2, xy=3' '\n']);
    if j<3
     Fext(i,j)=input(['Magnitud de la carga en el nudo ' num2str(i) '\n']);
    else 
     F=input(['Magnitud de la carga en el nudo ' num2str(i) '\n']);
     Ang=input(['Ángulo de inclinación de la carga con respecto a x positivo ' num2str(i) '\n']);
     Fext(i,1)=F*cosd(Ang);
     Fext(i,2)=F*sind(Ang);
    end   
end

%Secciones y materiales
for i=1:1:nsec;
    Secc(i,1)=input(['Área de la sección? ' num2str(i) '\n']);
    Secc(i,2)=input(['Inercia de la sección? ' num2str(i) '\n']);
end

%Modulo de elasticidad
for i=1:1:nmat;
    Mat(i,1)=input(['Módulo de elasticidad del material? ' num2str(i) '\n']);
end

%Identificación de los elementos
for i=1:1:nel
    Elem(i,1)=input(['Nudo inicial del elemento ' num2str(i) '\n']);
    Elem(i,2)=input(['Nudo final del elemento ' num2str(i) '\n']);
    Elem(i,3)=input(['Tipo de sección del elemento ' num2str(i) '\n']);
    Elem(i,4)=input(['Tipo de material del elemento ' num2str(i) '\n']);
    xi=Coord(Elem(i,1),1);
    yi=Coord(Elem(i,1),2);
    xf=Coord(Elem(i,2),1);
    yf=Coord(Elem(i,2),2);
    Delx=xf-xi;
    Dely=yf-yi;
    Elem(i,5)=sqrt(Delx^2+Dely^2);
    disp (Elem);
end

%Matriz de Transformación
for i=1:1:nel
    xi=Coord(Elem(i,1),1);
    yi=Coord(Elem(i,1),2);
    xf=Coord(Elem(i,2),1);
    yf=Coord(Elem(i,2),2);
    Delx=xf-xi;
    Dely=yf-yi;
    Long=sqrt(Delx^2+Dely^2);
    Cs=Delx/Long;
    Sn=Dely/Long;
    T=[Cs -Sn 0 0 0 0;Sn Cs 0 0 0 0;0 0 1 0 0 0;0 0 0 Cs -Sn 0;0 0 0 Sn Cs 0;0 0 0 0 0 1];
    Matrizdetransformacion(:,:,i) = T;
end

%Calculo de Fuerzas fijas
numerodeelemcaragdos = input ('cuantos elementos estan cargados? ');
for l = 1:1:numerodeelemcaragdos;
    inclelem = input('Elemento inclinado=1 o elemento ortogonal=2? ');
    if inclelem == 1;
        queelemento = input ('Cual elemento es el cargado? ');
        cargasinternas1 = input ('Ingrese el numero de cargas al interior del elemento inclinado');
        
        
        %for l=1:1:6
        %FF(GLel(queelemento,l),1)=FF(GLel(queelemento,l),1)+vfglob(l,1);        
        %end
    
        
        for i = 1:1:cargasinternas1;            
                Tipodecarga = input ('Transv distribuida=[1], Transv puntual=[2], Axial distribuida=[3], Axial puntual=[4] ');
                L = input('¿ Magnitud de la proyección? ');                
                if Tipodecarga == 1;
                w = input('¿ Magnitud y sentido de la carga transversal distribuida (kN/m)? ');
                vfglob(:,1,queelemento)=[0;-(w*L)/2;-(w*L^2)/12;0;-(w*L)/2;(w*L^2)/12];
                vfloc(:,1,queelemento)= Matrizdetransformacion(:,:,queelemento)*vfglob(:,1,queelemento);
                elseif Tipodecarga == 2;
                Puntual = input ('¿Magnitud y sentido de la carga transversal puntual (kN)? ');
                Distanciaizq = input ('¿Distancia de izquierda a derecha de la carga? ');
                Distanciader = input ('¿Distancia derecha izquierda de la carga? ');
                vfglob(:,1,queelemento)=[0;-((3*Distanciaizq)+(Distanciader))*((Puntual*Distanciader^2)/L^3);
                    -(Puntual*Distanciaizq*Distanciader^2)/L^2;0;-((Distanciaizq)+(3*Distanciader))*((Puntual*Distanciaizq^2)/L^3);
                    (Puntual*Distanciaizq^2*Distanciader)/L^2];
                vfloc(:,1,queelemento)= Matrizdetransformacion(:,:,queelemento)*vfglob(:,1,queelemento);
                elseif Tipodecarga == 3;
                w = input('¿ Magnitud y sentido de la carga axial distribuida (kN/m)? ');
                vfglob(:,1,queelemento)=[-(w*L)/2;0;0;-(w*L)/2;0;0];
                vfloc(:,1,queelemento)= Matrizdetransformacion(:,:,queelemento)*vfglob(:,1,queelemento);
                elseif Tipodecarga == 4;
                Puntual = input ('¿Magnitud y sentido de la carga axial puntual (kN)? ');
                Distanciaizq = input ('¿Distancia de izquierda a derecha de la carga? ');
                Distanciader = input ('¿Distancia derecha izquierda de la carga? ');
                vfglob(:,1,queelemento)=[-(Puntual*Distanciader)/L;0;0;-(Puntual*Distanciaizq)/L;0;0];
                vfloc(:,1,queelemento)= Matrizdetransformacion(:,:,queelemento)*vfglob(:,1,queelemento);
                end
        vf(:,1,queelemento)=vfglob(:,1,queelemento)+vf(:,1,queelemento);
        
        disp (vf);
        end
       
                
    elseif inclelem == 2;
                    queelemento = input ('Cual elemento es el cargado? ');
                    cargasinternas2 = input ('Ingrese el numero de cargas al interior del elemento ortogonal ');
                    for i = 1:1:cargasinternas2;
                        Tipodecarga = input ('Transv distribuida=[1], Transv puntual=[2], Axial distribuida=[3], Axial puntual=[4] ');
                        if Tipodecarga == 1;
                             w = input('¿ Magnitud y sentido de la carga transversal distribuida (kN/m)? ');
                            vfglob(:,1,queelemento)=[0;-(w*Elem(queelemento,5))/2;-(w*Elem(queelemento,5)^2)/12;
                                0;-(w*Elem(queelemento,5))/2;(w*Elem(queelemento,5)^2)/12];
                            vfloc(:,1,queelemento)= Matrizdetransformacion(:,:,queelemento)*vfglob(:,1,queelemento);
                        elseif Tipodecarga == 2;
                            Puntual = input ('¿Magnitud y sentido de la carga transversal puntual (kN)? ');
                            Distanciaizq = input ('¿Distancia de izquierda a derecha de la carga? ');
                            Distanciader = input ('¿Distancia derecha izquierda de la carga? ');
                            vfglob(:,1,queelemento)=[0;-((3*Distanciaizq)+(Distanciader))*((Puntual*Distanciader^2)/Elem(queelemento,5)^3);-(Puntual*Distanciaizq*Distanciader^2)/Elem(queelemento,5)^2;
                            0;-((Distanciaizq)+(3*Distanciader))*((Puntual*Distanciaizq^2)/Elem(queelemento,5)^3);(Puntual*Distanciaizq^2*Distanciader)/Elem(queelemento,5)^2];
                            vfloc(:,1,queelemento)= Matrizdetransformacion(:,:,queelemento)*vfglob(:,1,queelemento);
                        elseif Tipodecarga == 3;
                            w = input('¿ Magnitud y sentido de la carga axial distribuida (kN/m)? ');
                            vfglob(:,1,queelemento)=[-(w*Elem(queelemento,5))/2;0;0;-(w*Elem(queelemento,5))/2;0;0];
                            vfloc(:,1,queelemento)= Matrizdetransformacion(:,:,queelemento)*vfglob(:,1,queelemento);
                        elseif Tipodecarga == 4;
                            Puntual = input ('¿Magnitud y sentido de la carga axial puntual (kN)? ');
                            Distanciaizq = input ('¿Distancia de izquierda a derecha de la carga? ');
                            Distanciader = input ('¿Distancia derecha izquierda de la carga? ');
                            vfglob(:,1,queelemento)=[-(Puntual*Distanciader)/Elem(queelemento,5);0;0;-(Puntual*Distanciaizq)/Elem(queelemento,5);0;0];
                            vfloc(:,1,queelemento)= Matrizdetransformacion(:,:,queelemento)*vfglob(:,1,queelemento);                            
                        end
                    end
                        F(1) = MGL(Elem(queelemento,1),1);
                        F(2) = MGL(Elem(queelemento,2),1);
                        F(3) = MGL(Elem(queelemento,1),2);
                        F(4) = MGL(Elem(queelemento,2),2);
                        F(5) = MGL(Elem(queelemento,1),3);
                        F(6) = MGL(Elem(queelemento,2),3);
                        %for k=1:1:6;
                            %Femporamiento(F(k)) = Femporamiento(F(k))+ M(k);
                        %end
                        %disp (Femporamiento);  
    end
end
             
%Determinación de la matriz de rigidez y de transformación de cada elemento
for i=1:1:nel
    xi=Coord(Elem(i,1),1);
    yi=Coord(Elem(i,1),2);
    xf=Coord(Elem(i,2),1);
    yf=Coord(Elem(i,2),2);
    Delx=xf-xi;
    Dely=yf-yi;
    Long=sqrt(Delx^2+Dely^2);
    Cs=Delx/Long;
    Sn=Dely/Long;
    T=[Cs -Sn 0 0 0 0;Sn Cs 0 0 0 0;0 0 1 0 0 0;0 0 0 Cs -Sn 0;0 0 0 Sn Cs 0;0 0 0 0 0 1];
    Matrizdetransformacion(:,:,i) = T;
    EI=Mat(Elem(i,4),1)*Secc(Elem(i,3),2);
    EA=Mat(Elem(i,4),1)*Secc(Elem(i,3),1);
    %Matriz de rigidez local
    r11=EA/(Long);
    r22=12*EI/(Long^3);
    r23=6*EI/(Long^2);
    r33=4*EI/Long;    
    r36=2*EI/Long;
    kel=[r11 0 0 -r11 0 0;0 r22 r23 0 -r22 r23;0 r23 r33 0 -r23 r36;-r11 0 0 r11 0 0;
        0 -r22 -r23 0 r22 -r23;0 r23 r36 0 -r23 r33];
    kelocal(:,:,i)=kel;
    %Matriz de rigidez global
    keg=T*kel*T';
    kelglobal(:,:,i)=keg;
    %Identificación de grados de libertad por elemento
    for j=1:1:3
    GLel(i,j)=MGL(Elem(i,1),j);
    GLel(i,j+3)=MGL(Elem(i,2),j);
    end
    %Ensamblaje de la matriz de rigidez
    for l=1:1:6
        for m=1:1:6
            KT(GLel(i,l),GLel(i,m))=KT(GLel(i,l),GLel(i,m))+keg(l,m);
        end
    end
end

f = reshape(CAp',1,GL);
u = reshape(Fext',1,GL);
cont =0;
for i=1:1:GL
        if f(1,i) == 0    
            cont = 1+cont;
            Mfuezext(cont,1) = u(1,i);
        end
end

Desplazamiento = inv(KT(1:GIC,1:GIC))* Mfuezext;
Fuerzasdescon = KT(GIC+1:GL,1:GIC)*Desplazamiento;
Desplaztotal = zeros(GL,1);

for i = 1:1:GIC
    Desplaztotal(i,1)=Desplazamiento(i,1)+Desplaztotal(i,1);
end

%Construccion matrices de desplazamniento para cada elemento
for i=1:1:nel
    for j=1:1:6
        Vdesplaglobal(j,1,i)=Desplaztotal(GLel(i,j),1);
    end
end

%Construcción matrices de fuerzas globales para cada elemento
for i=1:1:nel
    FGlobelem(:,1,i)= kelglobal(:,:,i)*Vdesplaglobal(:,1,i);      
end

%Construcción matrices de fuerzas locales para cada elemento
for i=1:1:nel
    FLocelem(:,1,i)= Matrizdetransformacion(:,:,i)'*FGlobelem(:,1,i);      
end


