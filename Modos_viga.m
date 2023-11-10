function [w,x,U]=Modos_viga(n,bctype,bmpar,npoints)

% Programa para calcular Las frecuencias naturales y la forma del modo normalizado de masas para una 
% viga Euler-Bernoulli con una condición de frontera elegida.
% [w, x, U] = Modos_viga (n, bctype, bmpar, npoints) devolverá el enésimo
% modo de vibración a la frecuencia (w) y forma modal (U) de una viga de Euler-Bernoulli.
% Si n es un vector, devuelve las formas de modo correspondiente y
% frecuencias naturales
% Sin argumentos de salida, se dibujan los modos.
% Si solo se solicita un modo y no hay argumentos de salida, la forma  modal es animada.
% La condición de frontera se define de la siguiente manera:
% bctype = 1 Libre-Libre
% bctype = 2 Empotrada-Libre
% bctype = 3 Empotrada-Apoyo simple
% bctype = 4 Empotrada-Apoyo deslizante
% bctype = 5 Empotrado-Empotrado
% bctype = 6 Apoyo simple-Apoyo simple
% 
% Los parámetros de la viga se ingresan a través del vector bmpar:
% bmpar = [E I rho A L];
% 
% Ejemplo: viga de aluminio de 20 cm de largo con h = 1.5 cm, b = 3 cm
% para animar el cuarto modo para condiciones de contorno libre-libre
% correr los siguientes comandos
%
 %E=6e9;
 %I=(1/12)*0.03*0.005^3;
 %rho=2700;
 %A=0.03*0.005;
 %L=0.20;
 %Modos_viga(2,1,[E I rho A L]);

% %%%%%%%%%%%%%%%%%%%
 %bmpar = [E I rho A L];
% % Copyright Joseph C. Slater, 2007
% % Engineering Vibration Toolbox
% 

npoints=100;
E=bmpar(1);
I=bmpar(2);
rho=bmpar(3);
A=bmpar(4);
L=bmpar(5);

len=[0:(1/(npoints-1)):1]';  %Longitud normalizada de la viga

%Determinación de las frecuencias naturales y modos de vibración
%dependiendo de las condiciones de frontera dadas

if bctype==1
    desc=['Libre-Libre '];
    Bnllow=[0  0 4.73004074486 7.8532046241 10.995607838 14.1371654913 17.2787596574];
    for i=1:length(n)
        if n(i)>7
            %for i=6:n
            Bnl(i)=(2*n(i)-3)*pi/2;
            %end
        else
            Bnl(i)=Bnllow(n(i));

        end
    end
    for i=1:length(n)
        if n(i)==1
            w(i,1)=0;
            U(:,i)=1+len*0;
        elseif n(i)==2
            w(i,1)=0;
            U(:,i)=len-.5;
        else
            sig=(cosh(Bnl(i))-cos(Bnl(i)))/(sinh(Bnl(i))-sin(Bnl(i)));
            w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
            b=Bnl(i)*len;
            U(:,i)=cosh(b)+cos(b)-sig*(sinh(b)+sin(b));
        end
        
        %U(:,i)=U(:,i)/U(101,i);
    end
    

elseif bctype==2
    desc=['Empotrada-Libre '];
    Bnllow=[1.88 4.69 7.85 10.99 14.14];
    for i=1:length(n)
        if n(i)>5
            %for i=6:n
            Bnl(i)=(2*n(i)-1)*pi/2;
            %end
        else
            Bnl(i)=Bnllow(n(i));
        end
    end
    for i=1:length(n)
        sig=(sinh(Bnl(i))-sin(Bnl(i)))/(cosh(Bnl(i))-cos(Bnl(i)));
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        b=Bnl(i)*len;
        U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
        %U(:,i)=U(:,i)/U(101,i);
    end
    
elseif bctype==3
    desc=['Empotrada-Apoyo simple '];
    Bnllow=[3.93 7.07 10.21 13.35 16.49];
    for i=1:length(n)
        if n(i)>5
            %for i=6:n
            %Bnl(i)=(2*n(i)-1)*pi/2
            Bnl(i)=(4*n(i)+1)*pi/4;
            %end
        else
            Bnl(i)=Bnllow(n(i));
        end
    end

    for i=1:length(n)
        sig=(cosh(Bnl(i))-cos(Bnl(i)))/(sinh(Bnl(i))-sin(Bnl(i)));
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        b=Bnl(i)*len;
        U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
        %U(:,i)=U(:,i)/U(52,i);
    end

elseif bctype==4
    desc=['Empotrada-Deslizante '];
    Bnllow=[2.37 5.50 8.64 11.78 14.92];

    for i=1:length(n)
        if n(i)>5
            %for i=6:n
            %Bnl(i)=(2*n(i)-1)*pi/2
            Bnl(i)=(4*n(i)-1)*pi/4;
        else
            Bnl(i)=Bnllow(n(i));
        end
    end
    for i=1:length(n)
        sig=(sinh(Bnl(i))+sin(Bnl(i)))/(cosh(Bnl(i))-cos(Bnl(i)));
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        b=Bnl(i)*len;
        U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
        %U(:,i)=U(:,i)/U(101,i);
    end

elseif bctype==5
    desc=['Empotrada-Empotrada']
    Bnllow=[4.73 7.85 11 14.14 17.28];
    
        for i=1:length(n)
        if n(i)>5
            %for i=6:n
            %Bnl(i)=(2*n(i)-1)*pi/2
            Bnl(i)=(2*n(i)+1)*pi/2;
        else
            Bnl(i)=Bnllow(n(i));
        end
    end

    for i=1:length(n)
        sig=(cosh(Bnl(i))-cos(Bnl(i)))/(sinh(Bnl(i))-sin(Bnl(i)));
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        b=Bnl(i)*len;
        U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
        %U(:,i)=U(:,i)/U(52,i);
    end
elseif bctype==6
    desc=['Apoyo siple-Apoyo simple'];
    for i=1:length(n)
        Bnl(i)=n(i)*pi;
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        U(:,i)=sin(Bnl(i)*len);
    end
end

for i=1:length(n)
    U(:,i)=U(:,i)/sqrt(U(:,i)'*U(:,i)*rho*A*L);   
end


%stopstop=0;pausepause=0
global stopstop ppause;
ppause=0;
x=len*L;
%Rutina de graficación, si se requiere.
if nargout==0
    if length(n)~=1
        for i=1:length(n)
            plot(x,U(:,i))
            axis([0 L min(min(U)) max(max(U))])
            figure(gcf)
            title([desc,'  ','Modo ',int2str(i),'     Frecuencia natural = ',num2str(w(i)),' rad/s'])
            ylabel('Amplitud modal')
            xlabel('Longitud a lo largo de la barra - x')
            grid on
            disp('Presione enter para continuar')
            pause(0.1)
        end
    else
        nsteps=50;
        clf
        step=2*pi/(nsteps);
        i=0:step:(2*pi-step);
        hold off
        handle=uicontrol('style','pushbutton','units','normal','backgroundcolor','red','position', ...
            [0.94 .94 .05 .05],'String','Alto','callback','global stopstop;stopstop=1;');
        handle2=uicontrol('style','pushbutton','units','normal','backgroundcolor','yellow','position', ...
            [0.94 .87 .05 .05],'String','Pausa','callback','global ppause;ppause=1;');
        handle3=uicontrol('style','pushbutton','units','normal','backgroundcolor','green','position', ...
            [0.94 .80 .05 .05],'String','Reanudar','callback','global ppause;ppause=0;');
        
        stopstop=0;
        bb=0;
        while stopstop==0&bb<100
            bb=bb+1;
            for ii=[i  ]
                while ppause==1
                    pause(.01)
                    if stopstop==1
                        delete(handle), delete(handle2), delete(handle3)
                        return
                    end
                    
                end
                
                plot(x,U(:,1)*cos(ii))
                axis([0 L -max(abs(U)) max(abs(U))])
                grid on
                figure(gcf)
                title([desc,'  ','Modo ',int2str(n),'     \omega_n = ',num2str(w(1)/(2*pi)),' Hz'])
                ylabel('Amplitud modal')
                xlabel('Longitud a lo largo de la barra - x')
                drawnow
                %pause
            end
        end
        clear stopstop
        delete(handle), delete(handle2), delete(handle3)
    end
end

%Automatically check for updates
vtbchk
