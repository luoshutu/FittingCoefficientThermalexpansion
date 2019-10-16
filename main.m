
% ��ƷΪ��Ƭ����ο�����ɢ��ΪD = 85mm2/s������������L = 2.04mm
% FΪ����Ƶ��(Hz)��Ph_trueΪʵ�ʲ����λ��ֵ,Vol_trueΪʵ�ʲ���źŷ�ֵ(mV)
 L           =	2.04;
 F           =	[ 41     , 46    , 51    , 56    , 61    , 67    , 71    , 76    , 80    , 86    , 96    , 106   , 117   , 126   , 137   , 146   , 156   , 166   , 176   , 186   , 196   , 206   , 216   , 226   , 236   , 246   , 256   , 266   , 275   , 286   ];
 Vol_true    =	[ 1.98   , 1.77  , 1.62  , 1.53  , 1.40  , 1.31  , 1.22  , 1.16  , 1.10  , 1.04  , 0.95  , 0.85  , 0.79  , 0.73  , 0.67  , 0.61  , 0.61  , 0.55  , 0.52  , 0.49  , 0.40  , 0.49  , 0.46  , 0.43  , 0.43  , 0.40  , 0.40  , 0.37  , 0.37  , 0.34  ];
 Phase_true  =	-[23.14  , 28.70 , 33.37 , 37.00 , 40.24 , 43.45 , 45.55 , 47.89 , 49.87 , 51.48 , 54.35 , 56.01 , 58.63 , 59.48 , 60.95 , 62.79 , 62.74 , 64.18 , 64.23 , 64.29 , 65.24 , 66.80 , 67.07 , 68.20 , 67.38 , 68.63 , 67.75 , 70.02 , 69.15 , 70.71 ];

%% ��Ƶ��
%% ���ö��ַ��ķ�ʽ���
D           =   [1 (1+110)/2 110];
D_theory    =   0;
Phase       =   zeros(3,length(F));
numerator   =   zeros(3,length(F));
denominator =   zeros(3,length(F));
Error       =   zeros(3);
num         =   0;

while((D(3) - D(1)) > 0.1)
   for i = 1:3
       numerator(i,:)    =   1-((3*(sinh(L*sqrt((pi*F)/D(i))) + sin(L*sqrt((pi*F)/D(i)))))./((2*L*sqrt((pi*F)/D(i)).*(cosh(L*sqrt((pi*F)/D(i))) + cos(L*sqrt((pi*F)/D(i)))))));
       denominator(i,:)  =   (3*(sinh(L*sqrt((pi*F)/D(i))) - sin(L*sqrt((pi*F)/D(i)))))./((2*L*sqrt((pi*F)/D(i)).*(cosh(L*sqrt((pi*F)/D(i))) + cos(L*sqrt((pi*F)/D(i))))));
       Phase(i,:)        =   atan(-(numerator(i,:)./denominator(i,:)));
       Phase(i,:)        =   Phase(i,:).*(180/pi);      %����λ������ֵ��λ�ɻ���ת��Ϊ����
       Error(i)          =   sum((Phase(i,:) - Phase_true).^2)/length(F);  %��������
   end
   [~,num]    =   min(min(Error));
   D_theory   =   D(num);
   if(num == 1)                         %ȷ�����ַ���Χ
       D(3) = D(2);
       D(2) = (D(1) + D(3))/2;
   elseif(num == 2)
       if(Error(1) > Error(3))
           D(1) = D(2);
       else
           D(3) = D(2);
       end
       D(2) = (D(1) + D(3))/2;
   else
       D(1) = D(2);
       D(2) = (D(1) + D(3))/2;
   end
end

%% ��ͼ                                      
figure; 
plot( F, Phase(num,:), 'r' );   
hold on;
scatter( F, Phase_true, 's',  'b' );
set(gca,'FontSize',20);
legend('����ֵ','ʵ��ֵ');  title('��λ-Ƶ��ͼ');   xlabel('����Ƶ��(hz)');  ylabel('��λ(degree)');                      	
text(130,-50,['Dexp = ',num2str(D_theory),' mm2/s'],'FontSize',20);       %��ʾ����ɢ�ʵ���Ͻ��
text(70,-65,'��','FontSize',20);
axis tight; grid on;


%% ��Ƶ��
% �����׵�pΪ��Ч��Ȳ������˴�ȡֵΪ1��alphaΪ������Ʒ������������ϵ�������ҵ�����20���϶�������Ϊ23.2E-6/���϶�
% kΪ�����ȵ��ʣ��ο�ֵΪ2.35W/cm���϶ȡ�
p = 1; alpha = 2.32; k = 0.235;     %��λ����֮��

%% ���ö��ַ��ķ�ʽ���
D           =   [1 (1+110)/2 110];
D_theory    =   0;
Voltage     =   zeros(3,length(F));
numerator   =   zeros(3,length(F));
denominator =   zeros(3,length(F));
Error       =   zeros(3);
num         =   0;

while((D(3) - D(1)) > 0.1)
   for i = 1:3
       numerator(i,:)    =   1-((3*(sinh(L*sqrt((pi*F)/D(i))) + sin(L*sqrt((pi*F)/D(i)))))./((2*L*sqrt((pi*F)/D(i)).*(cosh(L*sqrt((pi*F)/D(i))) + cos(L*sqrt((pi*F)/D(i)))))));
       denominator(i,:)  =   (3*(sinh(L*sqrt((pi*F)/D(i))) - sin(L*sqrt((pi*F)/D(i)))))./((2*L*sqrt((pi*F)/D(i)).*(cosh(L*sqrt((pi*F)/D(i))) + cos(L*sqrt((pi*F)/D(i))))));
       Voltage(i,:)      =   abs(((p*alpha)./(k*L*((pi*F)/D(i)))).*sqrt(numerator(i,:).^2 + denominator(i,:).^2));
       Error(i)          =   sum((Voltage(i,:) - Vol_true).^2)/length(F);  %��������
   end
   [~,num]    =   min(min(Error));
   D_theory   =   D(num);
   if(num == 1)                         %ȷ�����ַ���Χ
       D(3) = D(2);
       D(2) = (D(1) + D(3))/2;
   elseif(num == 2)
       if(Error(1) > Error(3))
           D(1) = D(2);
       else
           D(3) = D(2);
       end
       D(2) = (D(1) + D(3))/2;
   else
       D(1) = D(2);
       D(2) = (D(1) + D(3))/2;
   end
end

%% ��ͼ                                      
figure; 
plot( F, F.*Voltage(num,:), 'r' );   
hold on;
scatter( F, F.*Vol_true, 's',  'b' );
set(gca,'FontSize',20);
legend('����ֵ','ʵ��ֵ');  title('����*Ƶ��-Ƶ��ͼ');   xlabel('����Ƶ��(hz)');  ylabel('����*Ƶ��(mv*hz)');                      	
text(150,85,['Dexp = ',num2str(D_theory),' mm2/s'],'FontSize',20);       %��ʾ����ɢ�ʵ���Ͻ��
text(150,84,'��','FontSize',20);
axis tight; grid on;
