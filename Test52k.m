close all
clear all
% close all
load Afinal
global FZ0 R0 KY

% global PCY1 PDY1 PDY2 PDY3 ...
%     PEY1 PEY2 PEY3 PEY4  ...
%     PKY1 PKY2 PKY3 ...
%     PHY1 PHY2 PHY3 ...
%     PVY1 PVY2 PVY3 PVY4;

global LFZO LCX LMUX LEX LKX  LHX LVX LCY LMUY LEY LKY LHY LVY ...
    LGAY LTR LRES LGAZ LXAL LYKA LVYKA LS LSGKP LSGAL LGYR A_final


% [name path] = uigetfile;
% addpath path
% load(name);
[FileName,PathName] = uigetfile(['E:\Mine\Tires\TTC_GUI''\*.mat'],'Select a tire file:');
input_filename = [PathName FileName];
load(input_filename)
pinterest = 12*6.89476;
SAinterest = 0;
vinterest = 40;
x = (abs(P-pinterest)<5) .* (abs(SA-SAinterest)<0.5) .* (abs(V-vinterest) <=3);
[inx iny] = find(x==1);
SA1 = SA(inx);
IA1 = IA(inx);
SL1 = SL(inx);
FY1 = FY(inx);
FX1 = FX(inx);
FZ1 = FZ(inx);
MZ1 = MZ(inx);
MX1 = MX(inx);
P1 = P(inx);
ET1 = ET(inx);
V1 = V(inx);
SL1 = SL(inx);
FX1 = FX(inx);
RL1 = RL(inx);


m=1:length(SL1);     % point counter
sp=csaps(m,SL1,.1); % fit a generic spline to locate zeros.

z=fnzeros(sp);      % location of zero crossings
z=round(z(1,:));   % no dups and integer indices

figure  % Just checkin'
plot(m,SL1,'r')
hold on
% plot(z,zeros(length(z)),'ko')
% xlim([0 3200])

z([2:2:length(z)])=[];  % drop kick the shutdown points;
plot(z,zeros(length(z)),'bo')
line([0 m(end)], [0 0],'color','k')
xlabel('Point Count')
ylabel('Slip Ratio')
legend('Test Data','Computed Slip Points of Interest'),legend Boxoff
clear fmdata
q=0;

for n = 1:(length(z)-1) % for some reason there are some repeat runs here. Skip them
    sa = SA1(z(n):z(n+1));
    fz = FZ1(z(n):z(n+1));
    fy = FY1(z(n):z(n+1));
    mz = MZ1(z(n):z(n+1));
    mx = MX1(z(n):z(n+1));
    rl = RL1(z(n):z(n+1));
    ia = IA1(z(n):z(n+1));
    fx = FX1(z(n):z(n+1));
    sl = SL1(z(n):z(n+1));

% fy= mean(fz)*ones(length(fz),1).*fy./fz;
% mz= mean(fz)*ones(length(fz),1).*mz./fz;
% mx= mean(fz)*ones(length(fz),1).*mx./fz;
% fz= mean(fz).*ones(length(fz),1);

%     [tmp,imn]=min(sa);
%     [tmp,imx]=max(sa);
%     p=1:length(sa);
%     rng=imx-50:imx+50;
%     pp=polyfit(p(rng),mz(rng)',3);
%     mzf=polyval(pp,p(rng));
% 
%     ind=find(abs(mzf-mz(rng)') > 7);
%     mz(rng(ind))=mzf(ind);
% 
%     rng=imn-50:imn+50;
%     pp=polyfit(p(rng),mz(rng)',3);
%     mzf=polyval(pp,p(rng));
% 
%     ind=find(abs(mzf-mz(rng)') > 7);
%     mz(rng(ind))=mzf(ind);
    
    sp_fyk = csaps(sl,fy,.1);
    sp_fx=csaps(sl,fx,.1);
    
    
    for temp=(min(sl)):0.01:(max(sl));
        q=q+1;
        fmdata(q,1)=temp;
        fmdata(q,2)=round(mean(ia));
        fmdata(q,3)=mean(fz);
        fmdata(q,4)=fnval(sp_fyk,temp);
        fmdata(q,5)=fnval(sp_fx,temp);
        fmdata(q,6)=mean(sa);
    end
end

fmdata  = sortrows(fmdata,[2,1,3]);

incls   = unique(round(fmdata(:,2)))'
nincls  = length(incls)

slips   =  unique(round(fmdata(:,1)))'
nslips  = length(slips)

t.SL= fmdata(:,1);
t.IA= fmdata(:,2);
t.FZ= fmdata(:,3);
t.FYK= fmdata(:,4);
t.FX= fmdata(:,5);
t.SA= fmdata(:,6);
%--------------------------------------------------------------------------
% Now patch in MF5.2 fitter:  

FZ0=  mean(t.FZ) % = FNOMIN = 'nominal wheel load'
R0  = .240

INPUT = [t.SA,t.FZ,t.IA,t.SL]; % slip, vert, incl

LFZO                     = 0.100000E+01         ;%typarr(  31)
LCX                      = 0.100000E+01         ;%typarr(  32)
LMUX                     = 0.100000E+01         ;%typarr(  33)
LEX                      = 0.100000E+01         ;%typarr(  34)
LKX                      = 0.100000E+01         ;%typarr(  35)
LHX                      = 0.100000E+01         ;%typarr(  36)
LVX                      = 0.100000E+01         ;%typarr(  37)
LCY                      = 0.100000E+01         ;%typarr(  38)
LMUY                     = 0.100000E+01         ;%typarr(  39)
LEY                      = 0.100000E+01         ;%typarr(  40)
LKY                      = 0.100000E+01         ;%typarr(  41)
LHY                      = 0.100000E+01         ;%typarr(  42)
LVY                      = 0.100000E+01         ;%typarr(  43)
LGAY                     = 0.100000E+01         ;%typarr(  44)
LTR                      = 0.100000E+01         ;%typarr(  45)
LRES                     = 0.100000E+01         ;%typarr(  46)
LGAZ                     = 0.100000E+01         ;%typarr(  47)
LXAL                     = 0.100000E+01         ;%typarr(  48)
LYKA                     = 0.100000E+01         ;%typarr(  49)
LVYKA                    = 0.100000E+01         ;%typarr(  50)
LS                       = 0.100000E+01         ;%typarr(  51)
LSGKP                    = 0.100000E+01         ;%typarr(  52)
LSGAL                    = 0.100000E+01         ;%typarr(  53)
LGYR                     = 0.100000E+01         ;%typarr(  54)

%[LATERAL_COEFFICIENTS]

% PCY1 =  -0.0088;
% PDY1 = -435.744;
% PDY2 = -8.181;
% PDY3 = 7.7211;
% PEY1 = 1.0626;
% PEY2 = -.0613;
% PEY3 = 0.0214;
% PEY4 = -0.2166;
% PKY1 = -53.2126;
% PKY2 = 2.1268;
% PKY3 = 2.1448;
% PHY1 = 0.0039;
% PHY2 =  0.0021;
% PHY3 = -.1602;
% PVY1 = -0.0146;
% PVY2 = -0.003;
% PVY3 = -3.3543;
% PVY4 = -1.1311;
RBY1				= +2.033e+001			;%typarr(109)
RBY2				= +8.152e+000			;%typarr(110)
RBY3				= -1.243e-002			;%typarr(111)
RCY1				= +9.317e-001			;%typarr(112)
REY1				= -3.982e-004			;%typarr(113)
REY2				= +3.077e-001			;%typarr(114)
RHY1				= +0.000e+000			;%typarr(115)
RHY2				= +0.000e+000			;%typarr(116)
RVY1				= +0.000e+000			;%typarr(117)
RVY2				= +0.000e+000			;%typarr(118)
RVY3				= +0.000e+000			;%typarr(119)
RVY4				= +0.000e+000			;%typarr(120)
RVY5				= +0.000e+000			;%typarr(121)
RVY6				= +0.000e+000			;%typarr(122)
clear AA RESNORM
A_str ={'RBY1' 'RBY2' 'RBY3' 'RCY1' 'REY1' 'REY2' 'RHY1' 'RHY2' 'RVY1' 'RVY2' 'RVY3' 'RVY4' 'RVY5' 'RVY6'};
A_old =[ RBY1   RBY2   RBY3   RCY1   REY1   REY2   RHY1   RHY2   RVY1   RVY2   RVY3   RVY4   RVY5   RVY6 ];

options =optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7);

fig1=figure ('MenuBar','none','Name',['Pacejka_97  fy Fitting Results'],'Position',[2 2 1600 1180],'NumberTitle','off');
for k=1:20
    [A,RESNORM(k),RESIDUAL,EXITFLAG] = lsqcurvefit('MF52_Fy_combined_fcn',A_old,INPUT,t.FYK,[],[],options);
        AA(:,k)=A; 
        for n=1:14
        subplot(3,6,n)
        bar([AA(n,:)],'group')
        title(['A(' num2str(n) ')' ' =' A_str{n}],'FontSize',8)
    end

    for n=1:14  % update A coefficients to newest values
        disp(['A_old(' num2str(n) ') = ' num2str(A_old(n)) ';    ' 'A(' num2str(n) ') = ' num2str(A(n)) ';'])
        eval(['A_old(' num2str(n) ') = ' num2str(A(n)) ' -1*eps*rand;']) % bootstrap
    end
%     set(fig1,'Name',[filename '     Free Rolling Lateral Force      Iteration: ' num2str(k)  '         RESNORM: ' num2str(RESNORM(k)) ])
    drawnow
end
pacefy=MF52_Fy_combined_fcn(A,INPUT);
inx0=find(INPUT(:,3)==0);
figure;
hold on