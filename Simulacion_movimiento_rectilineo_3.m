 clear 
 close all
 clc


t=0:0.001:6;
F1=exp(-2*t./(2*4.5)).*4.5.*cos(2*t);
F2=0.6*t.^2-0.1*t.^3;
f3=3.5-9.81*t.^2;
xt=F1;
vt=exp(-2*t./(2*4.5)).*(-4.5*2).*sin(2*t)+4.5*cos(2*t).*(-2/(2*4.5)).*exp(-2*t./(2*4.5));
at=exp(-2*t./(2*4.5)).*(-4.5*2*2).*cos(2*t)+(-4.5*2)*sin(2*t).*(-2/(2*4.5)).*exp(-2*t./(2*4.5))+4.5*cos(2*t).*((2*2)/(2*2*4.5*4.5)).*exp(-2*t./(2*4.5))-4.5*2*sin(2*t).*(-2/(2*4.5)).*exp(-2*t./(2*4.5));
%figure
%grid


%axis([-1.25*max(abs(xt)) 1.25*max(abs(xt)),-1.25*max(abs(xt)) 1.25*max(abs(xt))])
%axis([-5 5,-5 5])


  handle=uicontrol('style','pushbutton','units','normal','backgroundcolor','red','position', ...
            [0.85 .94 .13 .05],'String','Alto','callback','global stopstop;stopstop=1;'); 
  handle2=uicontrol('style','pushbutton','units','normal','backgroundcolor','yellow','position', ...
            [0.85 .87 .13 .05],'String','Pausa','callback','global ppause;ppause=1;');
  handle3=uicontrol('style','pushbutton','units','normal','backgroundcolor','green','position', ...
            [0.85 .8 .13 .05],'String','Reanudar','callback','global ppause;ppause=0;');

            


global stopstop ppause;
ppause=0;
i=1;
stopstop =0;

  while( i<6000&&stopstop==0)
    
     cla
     
     
     
     subplot(2,2,1)
     hold on
      plot([-5:0.001:5],0*[-5:0.001:5],'k -','LineWidth',2)
      plot(t(1:i),xt(1:i),'b -','LineWidth',2)
      axis([0 6,-5 5])
      set(gca,'defaulttextinterpreter','latex')
  
   set(get(gca,'XLabel'),'String','t [s]',...
                    'FontName','Arial',...
                    'FontAngle','normal',...
                    'FontSize',20)
% 
 set(get(gca,'YLabel'),'String','x(t) [m]',...
                     'FontName','Arial',...
                     'FontAngle','normal',...
                     'FontSize',20)

 set(gca,'fontsize',20);
set(gca,'fontname','Arial','FontWeight','normal');  
title(strcat('t= ',num2str((i-1)*0.001),' [s]','    x(t)= ',num2str(xt(i)),' [m]'))
    
   
   subplot(2,2,2)
     hold on
     plot([-5:0.001:5],0*[-5:0.001:5],'k -','LineWidth',2)
     plot(t(1:i),vt(1:i),'b -','LineWidth',2)
     axis([0 6,-10 10])
   set(gca,'defaulttextinterpreter','latex')
  
   set(get(gca,'YLabel'),'String','v(t) [m/s]',...
                    'FontName','Arial',...
                    'FontAngle','normal',...
                    'FontSize',20)

   set(get(gca,'XLabel'),'String','t [s]',...
                    'FontName','Arial',...
                    'FontAngle','normal',...
                    'FontSize',20)
 
 
 
set(gca,'fontsize',20);
set(gca,'fontname','Arial','FontWeight','normal');  
                
  title(strcat('t= ',num2str((i-1)*0.001),' [s]','    v(t)= ',num2str(vt(i)),' [m/s]'))
                  
                   subplot(2,2,3)
     hold on
     plot([-5:0.001:5],0*[-5:0.001:5],'k -','LineWidth',2)
      plot(t(1:i),at(1:i),'b -','LineWidth',2)
     %rectangle('Position',[xt(i) f3(i) 0.18 0.18],'FaceColor',[0 .5 .5],'Curvature',0.5)
     axis([0 6,-50 50])
   set(gca,'defaulttextinterpreter','latex')
  
   set(get(gca,'XLabel'),'String','t [s]',...
                    'FontName','Arial',...
                    'FontAngle','normal',...
                    'FontSize',20)

        set(get(gca,'YLabel'),'String','a(t) [m/s^2]',...
                    'FontName','Arial',...
                    'FontAngle','normal',...
                    'FontSize',20)             
                    
 
set(gca,'fontsize',20);
set(gca,'fontname','Arial','FontWeight','normal');  
      title(strcat('t= ',num2str((i-1)*0.001),' [s]','    a(t)= ',num2str(at(i)),' [m/s^2]'))           
                
                                  subplot(2,2,4)
     hold on
     plot([-5:0.001:5],0*[-5:0.001:5],'k -','LineWidth',2)
     rectangle('Position',[xt(i) 0.0 0.35 0.24],'FaceColor',[0 .5 .5],'Curvature',0.5)
     axis([-5 5,-1 1])
   set(gca,'defaulttextinterpreter','latex')
  
   set(get(gca,'XLabel'),'String','x(t) [m]',...
                    'FontName','Arial',...
                    'FontAngle','normal',...
                    'FontSize',20)

 
set(gca,'fontsize',20);
set(gca,'fontname','Arial','FontWeight','normal');  
                
         title(strcat('t= ',num2str((i-1)*0.001),' [s]'))           
                       
                while ppause==1
                    pause(.1) 
                    cla
                end
        i+=100;       
         pause(0.001)
       %cla
                  
end     
    
    




