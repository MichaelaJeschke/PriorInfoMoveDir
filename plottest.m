%histogram alle drei zero first middle last in einem
q = {'zero'}
 m = {'First'}
vektorallfirst = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})'];
m = {'Middle'}
vektorallmiddle = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})'];
m = {'Last'}
vektoralllast = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})'];

subAx1 = subplot(1, 3, 1, polaraxes);
subAx2 = subplot(1, 3, 2, polaraxes);
subAx3 = subplot(1, 3, 3, polaraxes);
obj1 = CircHist(vektorallfirst, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true, 'parent', subAx1);
obj2 = CircHist(vektorallmiddle, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true, 'parent', subAx2);
obj3 = CircHist(vektoralllast, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true, 'parent', subAx3);
obj1.setRLim(200, 20);
obj2.setRLim(200, 20);
obj3.setRLim(200, 20);
obj1.polarAxs.ThetaZeroLocation = 'right';
obj2.polarAxs.ThetaZeroLocation = 'right';
obj3.polarAxs.ThetaZeroLocation = 'right';

obj1.colorBar = [0 .4 .8];  
obj2.colorBar = [0 .4 .8];  
obj3.colorBar = [0 .4 .8];  
obj1.colorAvgAng = [1, 0.5, 0];
obj2.colorAvgAng = [1, 0.5, 0];
obj3.colorAvgAng = [1, 0.5, 0];

set(obj1.barHReflected, 'color', [1 1 1] * 0.6);
set(obj2.barHReflected, 'color', [1 1 1] * 0.6);
set(obj3.barHReflected, 'color', [1 1 1] * 0.6);

obj1.setThetaLabel('First', 'bottomleft');
obj2.setThetaLabel('Middle', 'bottomleft');
obj3.setThetaLabel('Last', 'bottomleft');

obj1.scaleBarSide = 'right'; 
obj2.scaleBarSide = 'right';
obj3.scaleBarSide = 'right';

%jetzt: FIRST ALLE QUALITIES GETRENNT

m = {'First'}
q = {'zero'}
vektorallzero = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})'];
q = {'fifteen'}
vektorallfifteen = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})'];
q = {'twentyfive'}
vektoralltwentyfive = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})'];
q = {'thirtyfive'}
vektorallthirtyfive = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})'];
q = {'fifty'}
vektorallfifty = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})'];



subAx1 = subplot(2, 3, 1, polaraxes);
subAx2 = subplot(2, 3, 2, polaraxes);
subAx3 = subplot(2, 3, 3, polaraxes);
subAx4 = subplot(2, 3, 4, polaraxes);
subAx5 = subplot(2, 3, 5, polaraxes);
obj1 = CircHist(vektorallzero, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true, 'parent', subAx1);
obj2 = CircHist(vektorallfifteen, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true, 'parent', subAx2);
obj3 = CircHist(vektoralltwentyfive, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true, 'parent', subAx3);
obj4 = CircHist(vektorallthirtyfive, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true, 'parent', subAx4);
obj5 = CircHist(vektorallfifty, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true, 'parent', subAx5);
obj1.setRLim(180, 20);
obj2.setRLim(180, 20);
obj3.setRLim(180, 20);
obj4.setRLim(180, 20);
obj5.setRLim(180, 20);
obj1.polarAxs.ThetaZeroLocation = 'right';
obj2.polarAxs.ThetaZeroLocation = 'right';
obj3.polarAxs.ThetaZeroLocation = 'right';
obj4.polarAxs.ThetaZeroLocation = 'right';
obj5.polarAxs.ThetaZeroLocation = 'right';

%BESCHRIFTUNG
obj1.setThetaLabel('0% quality', 'bottomleft');
obj2.setThetaLabel('15% quality', 'bottomleft');
obj3.setThetaLabel('25% quality', 'bottomleft');
obj4.setThetaLabel('35% quality', 'bottomleft');
obj5.setThetaLabel('50% quality', 'bottomleft');

%ZOOMEN ; QUALITIY LABELS VERSCHIEBEN UND FONT DER THETA ACHSE und löschen von dem oben MUSS MANUELL
%PASSIEREN


obj1.colorBar = [0 .4 .8];  
obj2.colorBar = [0 .4 .8];  
obj3.colorBar = [0 .4 .8];
obj4.colorBar = [0 .4 .8];  
obj5.colorBar = [0 .4 .8];
obj1.colorAvgAng = [1, 0.5, 0];
obj2.colorAvgAng = [1, 0.5, 0];
obj3.colorAvgAng = [1, 0.5, 0];
obj4.colorAvgAng = [1, 0.5, 0];
obj5.colorAvgAng = [1, 0.5, 0];
obj1.scaleBarSide = 'right'; 
obj2.scaleBarSide = 'right';
obj3.scaleBarSide = 'right';
obj4.scaleBarSide = 'right';
obj5.scaleBarSide = 'right';

histogrami = CircHist(vektorall, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true);
histogrami.setRLim(150, 20);
set(histogrami.barHReflected, 'color', [1 1 1] * 0.6);
histogrami.polarAxs.ThetaZeroLocation = 'right';
%ACHTUNG HIER NUN NOCH A::E RADIALEN ACHSEN G:EICH MACHEN UND HSTOGRAMM
%DREHEN NACH RECHTS AM BESTEN UND GESPEIGELTE BESCHRIFTUNG.....

%toPng([histogrami])

%close all
%close all hidden  
 %       end
%end





% für plot und anova mache ich eine vektor mit prozentanteil orthogonal des
% first move und gucke mir dan die unterschiede zwischen den qualities an.
for  vp = 1:16
     for trial=1:72     
 if expdata.subj(vp).shifted.zero.First(trial)>=75 && expdata.subj(vp).shifted.zero.First(trial)<=105 
     trialOrthzero(vp,trial)=1;
 else trialOrthzero(vp,trial)=0;
 end
 
 if expdata.subj(vp).shifted.fifty.First(trial)>=75 && expdata.subj(vp).shifted.fifty.First(trial)<=105 
     trialOrthfifty(vp,trial)=1;
 else trialOrthfifty(vp,trial)=0;
 end
 
 
 if expdata.subj(vp).shifted.fifteen.First(trial)>=75 && expdata.subj(vp).shifted.fifteen.First(trial)<=105 
     trialOrthfifteen(vp,trial)=1;
 else trialOrthfifteen(vp,trial)=0;
 end
 
 
 if expdata.subj(vp).shifted.twentyfive.First(trial)>=75 && expdata.subj(vp).shifted.twentyfive.First(trial)<=105 
     trialOrthtwentyfive(vp,trial)=1;
 else trialOrthtwentyfive(vp,trial)=0;
 end
 
 
 if expdata.subj(vp).shifted.thirtyfive.First(trial)>=75 && expdata.subj(vp).shifted.thirtyfive.First(trial)<=105 
     trialOrththirtyfive(vp,trial)=1;
 else trialOrththirtyfive(vp,trial)=0;
 end
 
 end
   
end

for i = 1:16
prozentvektorfirstzero(i, 1) = mean(trialOrthzero(i,:));
prozentvektorfirstfifty(i, 1) = mean(trialOrthfifty(i,:));
prozentvektorfirstfifteen(i, 1) = mean(trialOrthfifteen(i,:));
prozentvektorfirsttwentyfive(i, 1) = mean(trialOrthtwentyfive(i,:));
prozentvektorfirstthirtyfive(i, 1) = mean(trialOrththirtyfive(i,:));

end


subjID2 = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]'
 t3 = table(subjID2, prozentvektorfirstzero, prozentvektorfirstfifteen, prozentvektorfirsttwentyfive, prozentvektorfirstthirtyfive, prozentvektorfirstfifty)


zerocolumn = repmat(0,16,1);

fifteencolumn = repmat(15,16,1);

twentyfivecolumn = repmat(25,16,1);

thirtyfivecolumn = repmat(35,16,1);

fiftycolumn = repmat(50,16,1);

prozentvektorfirstzero = [prozentvektorfirstzero zerocolumn];
prozentvektorfirstfifty = [prozentvektorfirstfifty  fiftycolumn];
prozentvektorfirstfifteen = [prozentvektorfirstfifteen fifteencolumn];
prozentvektorfirsttwentyfive = [prozentvektorfirsttwentyfive twentyfivecolumn];
prozentvektorfirstthirtyfive = [prozentvektorfirstthirtyfive thirtyfivecolumn];

allprozentvektor = vertcat(prozentvektorfirstzero, prozentvektorfirstfifteen, prozentvektorfirsttwentyfive, prozentvektorfirstthirtyfive, prozentvektorfirstfifty);
plot(allprozentvektor(:,2), allprozentvektor(:,1))



%HIER BEGIUNNT ANOVAs

qualityvektor = allprozentvektor(:,2);
allprozentvektor(:,2)= [];
subjID = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]';

%package anova ienfaktirell messwiedeholt. gute

qualityletzte = [1; 1; 1; 1; 1; 1; 1;1; 1; 1; 1; 1; 1; 1; 1;1; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 3; 3; 3; 3; 3; 3; 3; 3; 3; 3; 3; 3; 3; 3; 3; 3; 4; 4; 4; 4; 4; 4; 4; 4; 4; 4; 4; 4; 4; 4; 4; 4; 5;  5; 5; 5; 5; 5; 5; 5; 5; 5; 5; 5; 5; 5; 5; 5];

anovamatrixorthoalt= [allprozentvektor, qualityletzte, subjID]

%[RMAOV1] = RMAOV1(anovamatrixortho)

for i= 1:16
   strokevektor = [expdata.subj(i).numstrokes, expdata.subj(i).quality], 
   strokevektor = sortrows(strokevektor, 2);
   strokevektorzero(i,1)= mean(strokevektor(1:72));
   strokevektorfifteen(i,1) = mean(strokevektor(73:144));
   strokevektortwentyfive(i,1) = mean(strokevektor(145:216));
   strokevektorthirtyfive(i,1) = mean(strokevektor(217:288));
   strokevektorfifty(i,1) = mean(strokevektor(289:360));
   
end

%kurz drescripitves

vektorprozentortho = [prozentvektorfirstzero, prozentvektorfirstfifteen, prozentvektorfirsttwentyfive, prozentvektorfirstthirtyfive, prozentvektorfirstfifty];
vektorprozentortho(:,[2 4 6 8 10]) = []
for i =1:5
   
    orthooverview.means(i) =  mean(vektorprozentortho(:,i))
    orthooverview.sd(i) = std(vektorprozentortho(:,i))
end



%jetzt dasselkbe für die strokes 
%mean strokes 



strokevektorfinal = [strokevektorzero; strokevektorfifteen; strokevektortwentyfive;  strokevektorthirtyfive; strokevektorfifty];
   
   anovamatrixstroke= [strokevektorfinal, qualityletzte, subjID]
   
   %[RMAOV1] = RMAOV1(anovamatrixstroke)

   
   %nochmal plot WORKS LIKE THAT NUR LIMITS UND EVtL HJINKRIEGEN DASS
   %PUBNKTE DRIN ODER SO
   
   %für r legngth plot zweipaltige matrix mit links quality und rechts r
   %werten
   
   rlengthplot = zeros(80,2);
   rlengthplot(:,1) = qualityvektor;
   rlengthplot(:,2) = [zeroR; fifteenR; twentyfiveR; thirtyfiveR; fiftyR]
   
   g=gramm('x',rlengthplot(:,1),'y',rlengthplot(:,2),'color', rlengthplot(:,1));
%mit whatever ausprobieren
   g.stat_summary('geom',{'bar','black_errorbar'});
g.draw()

%ortho

orthoplot = [qualityv ektor, allprozentvektor]
   g=gramm('x',orthoplot(:,1),'y',orthoplot(:,2),'color',orthoplot(:,1));
    g=gramm('x',Results_VP1(:,2),'y',Results_VP1(:,3),'color', Results_VP1(:,2));
 %  g.set_names('x','Quality of PVI','y','Orthogonal strokes (in %)');
g.stat_boxplot();
g.set_title('stat_boxplot()');
 g.draw()
   