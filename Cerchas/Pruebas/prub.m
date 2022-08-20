%Datos iniciales de entrada
clear
clear all
clc

nnd=3;
nel=3;
nsec=1;
nmat=1;
GIC=3;
GL=nnd*2;
nrest=GL-GIC;
nFext=1;
%Dimensionamiento de matrices
Coord=zeros(nnd,2);
CAp=zeros(nnd,2);
Secc=zeros(nsec,1);
Mat=zeros(nmat,1);
Fext=zeros(nnd,2);
Elem=zeros(nel,4);
MGL=zeros(nnd,2);
KT=zeros(GL,GL);
GLel=zeros(nel,4);

%Entrada de datos de los nudos
for i=1:1:nnd;
    nudo=i;
    Coord(i,1)=input(['Coordenada x del nudo ' num2str(i) '\n']);
    Coord(i,2)=input(['Coordenada y del nudo ' num2str(i) '\n']);
end

%Entrada de datos de las restricciones
for ires=1:1:nrest;
    i=input(['Nudo restringido:' '\n']);
    j=input(['restringido en x=1, o en y=2?? ']);
    CAp(i,j)=1;
end
%Matriz de grados de libertad
cont1=1;
cont2=GIC+1;
for i=1:1:nnd
    for j=1:1:2
        if CAp(i,j)==1
            MGL(i,j)=cont2;
            cont2=cont2+1;
        else
            MGL(i,j)=cont1;
            cont1=cont1+1;
        end
    end
end
%Secciones y materiales
for i=1:1:nsec;
    Secc(i,1)=input(['�rea de la secci�n ' num2str(i) '\n']);
end
for i=1:1:nmat;
    Mat(i,1)=input(['M�dulo de elasticidad del material ' num2str(i) '\n']);
end
   
%Fuerzas externas
for cont1=1:1:nFext
    i=input(['Nudo con carga' '\n']);
    j=input(['Direcci�n de la carga en x=1, y=2, xy=3' '\n']);
    if j<3
     Fext(i,j)=input(['Magnitud de la carga en el nudo ' num2str(i) '\n']);
    else 
     F=input(['Magnitud de la carga en el nudo ' num2str(i) '\n']);
     Ang=input(['�ngulo de inclinaci�n de la carga con respecto a x positivo ' num2str(i) '\n']);
     Fext(i,1)=F*cosd(Ang);
     Fext(i,2)=F*sind(Ang);
    end   
end

%Identificaci�n de los elementos
for i=1:1:nel
    Elem(i,1)=input(['Nudo inicial del elemento ' num2str(i) '\n']);
    Elem(i,2)=input(['Nudo final del elemento ' num2str(i) '\n']);
    Elem(i,3)=input(['Tipo de secci�n del elemento ' num2str(i) '\n']);
    Elem(i,4)=input(['Tipo de material del elemento ' num2str(i) '\n']);
end

%Determinaci�n de la matriz de rigidez y de transformaci�n de cada elemento
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
    T=[Cs -Sn 0 0;Sn Cs 0 0; 0 0 Cs -Sn;0 0 Sn Cs];
    Matrizdetransformacion(:,:,i) = T;
    AE=Secc(Elem(i,3),1)*Mat(Elem(i,4),1);
    %Matriz de rigidez local
    kel=(AE/Long)*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0];
    kelocal(:,:,i)=kel;
    %Matriz de rigidez global
    keg=T*kel*T';
    kelglobal(:,:,i)=keg;
    %Identificaci�n de grados de libertad por elemento
    for j=1:1:2
    GLel(i,j)=MGL(Elem(i,1),j);
    GLel(i,j+2)=MGL(Elem(i,2),j);
    end
    %Ensamblaje de la matriz de rigidez
    
    for l=1:1:4
        for m=1:1:4
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
Fuerzasdescom = KT(GIC+1:GL,1:GIC)*Desplazamiento;

Desplaztotal = zeros(GL,1);

for i = 1:1:GIC
    Desplaztotal(i,1)=Desplazamiento(i,1)+Desplaztotal(i,1);
end


%Construccion matrizes de desplazamniento para cada elemento

for i=1:1:nel
    for j=1:1:4
        matricitasdesplazglobal(j,1,i)=Desplaztotal(GLel(i,j),1);
    end
end

%Construccion Fuerzas en cada elemento Globalmente

for i=1:1:nel
    Fuerzasglobaleselemtnos(:,:,i) = kelglobal(:,:,i)*matricitasdesplazglobal(:,:,i);
end

%Construccion Fuerzas en cada elemento localmente
for i=1:1:nel
    Fuerzaslocaleselementos(:,:,i) = Matrizdetransformacion(:,:,i)'*Fuerzasglobaleselemtnos(:,:,i);
end

    









