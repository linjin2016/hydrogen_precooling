clear
clc
%-----------------Initial design data------------------------

%Temperature flow pressure parameters of hot fluid, hydrogen
M_h=120/3600; %  total mass flow rate kg/s
T_h_in=40+273.15; % inlet temperature K
T_h_out=-30+273.15; %outlet temperature K
P_h_in=90e3; % inlet pressure kPa

%Temperature flow pressure parameters of cold fluid, LM-8
M_c=9000/3600; % total mass flow rate kg/s
T_c_in=-40+273.15; % inlet temperature K
P_c_in=1.5e3; % inlet pressure kPa

%node number
N=30; %

%----------------------Total heat transfer------------------
Hhin=refpropm('H','T',T_h_in,'P',P_h_in,'HYDROGEN');
Hhout=refpropm('H','T',T_h_out,'P',P_h_in,'HYDROGEN');
Qtot=M_h*(Hhin-Hhout);
qseg=Qtot/N;

%----------------------structure parameter------------------
%Hot side 
Dch=0.89e-03;% Diameter of hot side channel m
twh=1.27e-03;% Thickness of hot side plate m
trh=1.66e-03;%Hot side plate rib width m
Ah=(pi*(Dch/2)^2)/2; % The cross-sectional area of a single semicircular channel on the hot side m2
Deh=pi*Dch/2/(1+pi/2); % Hydraulic diameter of a single semicircular channel m

%Cold side
ef=5.98e-03;%Plate-fin channel height m
Pt=2.21e-03; %Plate-fin channel width m
delt=0.29e-03; %Plate fin plate-thickness m
tt=0.78e-03; %separator-plate thickness m
lss=9.6e-03; %Length of a single serrated fin m

alfack=Pt/ef;
belt=delt/lss;
gamma=delt/Pt;

Ac=Pt*ef; % The cross-sectional area of a single semicircular channel on the cold side m2
Dec=4*ef*Pt*lss/(2*(Pt*lss+ef*lss+delt*ef)+delt*Pt);%Hydraulic diameter of a single plate-fin channel

delttwh=(twh-pi/8*Deh)*1e-03;%Hot side equivalent plate thickness m
delttwc=tt;%Cold side equivalent plate thickness m
delttw=(delttwh+delttwc)/2;%Arithmetic average of equivalent thickness on both sides of cold and hot m

Nph=23; %Number of hot side plates
Npc=Nph; %The number of cold side panels is the same as that of hot side panels
Nch=21;%Number of channels in a single plate on the hot side
Wtot=Nch*(Dch+trh);%Total width of hot side plate
Ncc=Wtot/(Pt+delt);%Number of channels in a single plate on the cold side
Htot=twh*Nph+(ef+tt)*Npc;% Total height
conw=16.3;%Thermal conductivity of 316SS sheet material W/mK

m_h=M_h/Nph/Nch;% Mass flow of a single hot side channel
m_c=M_c/Npc/Ncc;% Mass flow of a single cold side channel

T_c_out=-35+273.15; %  Customize the initial value: outlet temperature of cold side K
P_c_out=1.5e3; % Customize the initial value: outlet pressure of cold side kPa

T_h=ones(1,N+1); 
T_c=ones(1,N+1);  
P_h=ones(1,N+1); 
P_c=ones(1,N+1);  
T_h_b=ones(1,N);  
T_c_b=ones(1,N);  
P_h_b=ones(1,N);  
P_c_b=ones(1,N);  
delt_P_h=ones(1,N);
delt_P_c=ones(1,N);
delt_T_h=ones(1,N);
delt_T_c=ones(1,N);

T_h(1)=T_h_in;
P_h(1)=P_h_in;
T_c(1)=T_c_out;
P_c(1)=P_c_out;

i=1;
m=1;
n=1;
a=1;
b=1;
c=1;
d=1;
errg=1;
errh=1;
errm=1;
errn=1;

 while abs(n-T_c_in)>0.001|abs(m-P_c_in)>0.0001 

    Sxuyong=117000;% 316 stainless steel allowable stress kPa
    JointF=0.7; %joint factor
    jione=Sxuyong*JointF;
    jitwo=1.5*Sxuyong*JointF;
    %---Hot side wall thickness check---% 
    Smthx2=P_h_in*Dch/2/(2*(twh-Dch/2));% membrane stress kPa
    Sbthx2=P_h_in*Dch*Dch*(twh-Dch/2)/2/((twh-Dch/2)*(twh-Dch/2)*(twh-Dch/2));%bending stress
    Stotthx2=Smthx2+Sbthx2;%total stress kPa

    if Smthx2>jione || Stotthx2>jitwo
    return;
    end
    
    Smths2=P_h_in*Dch/2/(2*(tt));% membrane stress kPa
    Sbths2=P_h_in*Dch*Dch*(tt)/2/((tt)*(tt)*(tt));%bending stress
    Stotths2=Smths2+Sbths2;%total kPa

    if Smths2>jione || Stotths2>jitwo
    return;
    end
    
    
    Smth3=P_h_in*Dch/trh;
    if Smth3>jione
    return;
    end
    
    %-------Cold side wall thickness check-----%
    Smtc2=P_c_in*ef/(2*delt);% membrane stress kPa
    Sbtc2=P_c_in*Pt*Pt*delt/(delt)*(delt)*(delt);%bending stress
    Stottc2=Smtc2+Sbtc2;%total stress kPa

    if Smtc2>jione || Stottc2>jitwo
    return;
    end

    Smtc3=P_c_in*Pt/delt;
    if Smtc3>jione
    return;
    end


    delt_P=1;
    delt_T=1;
    
for j=1:N
    delt_P_h(j)=delt_P;
    delt_P_c(j)=delt_P;
    delt_T_h(j)=delt_T; 
    delt_T_c(j)=delt_T; 
end
    
   i=1;
   
while i<N+1
    
    i=i+1;
   
    delt_T_h_new=delt_T_h(i-1);
    delt_T_c_new=delt_T_c(i-1);
    delt_P_h_new=delt_P_h(i-1);
    delt_P_c_new=delt_P_c(i-1);
    
    errg=delt_T_h_new;
    errh=delt_T_c_new;
    errm=delt_P_h_new;
    errn=delt_P_c_new;
    
    
    while abs(errg)>0.001|abs(errh)>0.001|abs(errm)>0.0001|abs(errn)>0.0001

    T_h(i)=T_h(i-1)-delt_T_h(i-1);
    T_c(i)=T_c(i-1)-delt_T_c(i-1);
    P_h(i)=P_h(i-1)-delt_P_h(i-1);
    P_c(i)=P_c(i-1)+delt_P_c(i-1);
    T_h_b(i-1)=(T_h(i)+T_h(i-1))/2;
    T_c_b(i-1)=(T_c(i)+T_c(i-1))/2;
    P_h_b(i-1)=(P_h(i)+P_h(i-1))/2;
    P_c_b(i-1)=(P_c(i)+P_c(i-1))/2;
    
    
    Cph(i-1)=refpropm('C','T',T_h_b(i-1),'P',P_h_b(i-1),'HYDROGEN'); 
    Miuh(i-1)=refpropm('V','T',T_h_b(i-1),'P',P_h_b(i-1),'HYDROGEN');
    Denh(i-1)=refpropm('D','T',T_h_b(i-1),'P',P_h_b(i-1),'HYDROGEN');
    Conh(i-1)=refpropm('L','T',T_h_b(i-1),'P',P_h_b(i-1),'HYDROGEN'); 
    Reh(i-1)=m_h*Deh/(Miuh(i-1)*Ah); 
    Prh(i-1)=refpropm('^','T',T_h_b(i-1),'P',P_h_b(i-1),'HYDROGEN'); 
  
    Cpc(i-1)=2574.3; 
    Miuc(i-1)=2.98e-03;
    Denc(i-1)=1383;
    Conc(i-1)=0.415;
    
    Rec(i-1)=m_c*Dec/(Miuc(i-1)*Ac); 
    Prc(i-1)=Cpc(i-1)*Miuc(i-1)/Conc(i-1);
    
    %Hot side flow condition judgment---%
    %---laminar flow
    if Reh(i-1)<=2300
        Nuh(i-1)=4.089;
        fh(i-1)=63.12/Reh(i-1);
    end
    %---Sufficient turbulent condition
    if Reh(i-1)>=3000&&Reh(i-1)<5e6
    fh(i-1)=1/((1.82*log10(Reh(i-1))-1.64).^2);% f
    Nuh(i-1)=fh(i-1)/8*(Reh(i-1)-1000)*Prh(i-1)/(1+12.7*sqrt(fh(i-1)/8)*(Prh(i-1)^(2/3)-1));% Nu
    end
    %--transition section
    if Reh(i-1)>2300&&Reh(i-1)<3000
     fhshou=63.12/2300;
     Nuhshou(i-1)=4.089;
     fhmo=(1.82*log10(3000)-1.64)^-2;
     Nuhmo(i-1)=(fhmo/8)*(3000-1000)*Prh(i-1)/(1+12.7*sqrt(fhmo/8)*(Prh(i-1)^(2/3)-1));
     Nuh(i-1)=(Nuhmo(i-1)-Nuhshou(i-1))*(Reh(i-1)-2300)/(3000-2300)+Nuhshou(i-1);
     fh(i-1)=(fhmo-fhshou)*(Reh(i-1)-2300)/(3000-2300)+fhshou;
    end   
    hh(i-1)=Nuh(i-1)*Conh(i-1)/Deh;
    
    %cold side flow condition judgment---%
    fcf(i-1)=9.6243*Rec(i-1)^(-0.7422)*alfack^(-0.1856)*belt^0.3053*gamma^(-0.2659)*(1+7.669*0.00000001*Rec(i-1)^4.429*alfack^0.920*belt^3.767*gamma^0.236)^0.1;
    jc(i-1)=0.6522*Rec(i-1)^(-0.5403)*alfack^(-0.1541)*belt^0.1499*gamma^(-0.0678)*(1+5.269*0.00001*Rec(i-1)^1.34*alfack^0.504*belt^0.456*gamma^(-1.055))^0.1;

    fc(i-1)=4*fcf(i-1);
    Nuc(i-1)=jc(i-1)*Rec(i-1)*Prc(i-1)^(1/3);
    
    hc(i-1)=Nuc(i-1)*Conc(i-1)/Dec;%                                                                                                计算冷侧的表面对流传热系数
    
    %Overall heat transfer coefficient
    Atoth=(Dch+pi*Dch/2)*Nph*Nch;%The total heat transfer area of the hot side of the micro segment m2
    Atotc=(ef*2+Pt*2)*Npc*Ncc;%The total heat transfer area of the cold side of the micro segment m2
    hall(i-1)=1/(1/hh(i-1)+Atoth/hc(i-1)/Atotc+delttw/conw*Atoth/Atotc);% The overall heat transfer coefficient based on the hot side area
    
    L(i-1)=qseg/(hall(i-1)*Atoth*(T_h_b(i-1)-T_c_b(i-1))); %The length of a local bit segment,m

    delt_T_h_new=qseg/Cph(i-1)/M_h;
    delt_T_c_new=qseg/Cpc(i-1)/M_c;
    delt_P_h_new=fh(i-1)*L(i-1)*m_h^2/(Denh(i-1)*Deh*(Ah)^2)/2/1000; %kPa
    delt_P_c_new=fc(i-1)*L(i-1)*m_c^2/(Denc(i-1)*Dec*(Ac)^2)/2/1000; %kPa
    
    errg=delt_T_h_new-delt_T_h(i-1);
    errh=delt_T_c_new-delt_T_c(i-1);
    errm=delt_P_h_new-delt_P_h(i-1);
    errn=delt_P_c_new-delt_P_c(i-1);
    
    if abs(errg)>0.001|abs(errh)>0.001
        delt_T_h(i-1)=delt_T_h_new;
        delt_T_c(i-1)=delt_T_c_new;
    end
   
    if abs(errm)>0.0001|abs(errn)>0.0001
        delt_P_h(i-1)=delt_P_h_new;
        delt_P_c(i-1)=delt_P_c_new;
    end 
      
    end
    
end
  
       n=T_c(N+1);
       m=P_c(N+1);
       
       if abs(n-T_c_in)>0.001
          c=T_c(N+1)-T_c_in;
          dq=abs(0.5*c);
          if c>0
              d=T_c(1)-dq;
          elseif c<0
              d=T_c(1)+dq;
          end
          T_c(1)=d;
     end
     
      if abs(m-P_c_in)>0.0001
          a=P_c(N+1)-P_c_in;
          dp=abs(0.5*a);
          if a>0
              b=P_c(1)-dp;
          elseif a<0
              b=P_c(1)+dp;
          end
          P_c(1)=b;          
    end

 end
 
 Ltot=sum(L);
 
 %total pressure
 Denhshou=refpropm('D','T',T_h_b(1),'P',P_h_b(1),'HYDROGEN');%Hot side inlet density
 Denhmo=refpropm('D','T',T_h_b(N),'P',P_h_b(N),'HYDROGEN');%Hot side outlet density
 sudushou=m_h/Denhshou/Ah;
 sudumo=m_h/Denhmo/Ah;
 Ph_loss_ac=(-Denhmo*sudumo*sudumo+Denhshou*sudushou*sudushou)/1000;%kPa;
 Ph_loss_f=P_h(1)-P_h(N+1);
 
 Pc_loss_tot=P_c(N+1)-P_c(1);
 Ph_loss_tot=Ph_loss_f+Ph_loss_ac;
 P_loss_tot=Pc_loss_tot+Ph_loss_tot;
 
 V=((Dch+trh)*Nch)*((twh+ef+delt+tt)*Nph)*Ltot;
 Kall=sum(hall)/N;

 %hot side
 Smthx2=P_h_in*Dch/2/(2*(twh-Dch/2));%  kPa
 Sbthx2=P_h_in*Dch*Dch*(twh-Dch/2)/2/((twh-Dch/2)*(twh-Dch/2)*(twh-Dch/2));%
 Stotthx2=Smthx2+Sbthx2;% kPa

 Smths2=P_h_in*Dch/2/(2*(tt));%  kPa
 Sbths2=P_h_in*Dch*Dch*(tt)/2/((tt)*(tt)*(tt));%
 Stotths2=Smths2+Sbths2;% kPa

 Smth3=P_h_in*Dch/trh;
 
 %cold side
 Smtc2=P_c_in*ef/(2*delt);%  kPa
 Sbtc2=P_c_in*Pt*Pt*delt/(delt)*(delt)*(delt);
 Stottc2=Smtc2+Sbtc2;% kPa
 Smtc3=P_c_in*Pt/delt;
 
 %synergistic angle between the local heat transfer coefficient and local heat transfer temperature difference
kallsum=sqrtm(sum(hall(:).^2));
deltTlengre=T_h_b-T_c_b;
deltTlengresum=sqrtm(sum(deltTlengre(:).^2));
thjuzhenji=hall.*deltTlengre;
thjuzhenhe=sum(thjuzhenji(:));
cosalfa=thjuzhenhe/kallsum/deltTlengresum;
alfa=acosd(cosalfa);

Recave=sum(Rec)/N;
Rehave=sum(Reh)/N;

Areachuanreh=(pi*Dch/2+Dch).*L;
Areachuanrec=(2*(ef+Pt)).*L;
rezuhc=1./(hc.*Areachuanrec);
rezuhh=1./(hh.*Areachuanreh);
rezuconh=sum(rezuhh); % thermal resistance
rezuconc=sum(rezuhc);% thermal resistance

Areachuanreh=(pi*Dch/2+Dch).*L;%
Areachuanrec=(2*(ef+Pt)).*L;
rezuhcfenbu=1./(hc.*Areachuanrec);
rezuhhfenbu=1./(hh.*Areachuanreh);
rezuhc=1/(mean(hc)*(2*(ef+Pt))*Ltot);
rezuhh=1/(mean(hh)*(pi*Dch/2+Dch)*Ltot);
rezuzong=1/(mean(hall)*(pi*Dch/2+Dch)*Ltot);

resulttwo=[Ltot,Pc_loss_tot,Ph_loss_tot,m_c,m_h,Recave,Rehave,V,alfa,rezuhc,rezuhh,rezuzong];
resultthree=[Dec,Ltot,Pc_loss_tot,Ph_loss_tot,m_c,m_h,Recave,Rehave,V,alfa,rezuhc,rezuhh,rezuzong]; 
resultli=[Smthx2/1000,Stotthx2/1000,Smth3/1000,Smtc2/1000,Stottc2/1000,Smtc3/1000,Smths2/1000,Stotths2/1000];