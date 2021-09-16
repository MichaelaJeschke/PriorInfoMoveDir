%erste erstes viertel zweite zweites dritte drittes vierte viertes und fünfte sind alle
%so für alle 15 dinger

   
move = {'First', 'Middle', 'Last'};    
qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'};    
p = [1:16];
vp = 0;
for person = p
    vp = vp+1;
    for q = qua
        for m = move
    testrad = circ_ang2rad(expdata.subj(vp).shifted.(q{1}).(m{1})');
     [r, mu] = circ_axialmean(testrad(1:17, 1), 2, 1); 
     means.(m{1}).(q{1}){vp,1} = mu;
     testrad = circ_ang2rad(expdata.subj(vp).shifted.(q{1}).(m{1})');
     [r, mu] = circ_axialmean(testrad(18:36, 1), 2, 1); 
     means.(m{1}).(q{1}){vp,2} = mu;
     testrad = circ_ang2rad(expdata.subj(vp).shifted.(q{1}).(m{1})');
     [r, mu] = circ_axialmean(testrad(37:55,1), 2, 1); 
     means.(m{1}).(q{1}){vp,3} = mu;
     testrad = circ_ang2rad(expdata.subj(vp).shifted.(q{1}).(m{1})');
     [r, mu] = circ_axialmean(testrad(56: 72,1), 2, 1); 
     means.(m{1}).(q{1}){vp,4} = mu;
     testrad = circ_ang2rad(expdata.subj(vp).shifted.(q{1}).(m{1})');
     [r, mu] = circ_axialmean(testrad, 2, 1); 
     means.(m{1}).(q{1}){vp,5} = mu;
        end
    end
end

%jetzt die negativen plus pi - noch loopen


for i = 1:16
  for q = qua
        for m = move
    if means.(m{1}).(q{1}){i,1} < 0
        means.(m{1}).(q{1}){i,1} = means.(m{1}).(q{1}){i,1} + pi;
    else means.(m{1}).(q{1}){i,1} = means.(m{1}).(q{1}){i,1};
    end
    
     if means.(m{1}).(q{1}){i,2} < 0
        means.(m{1}).(q{1}){i,2} =  means.(m{1}).(q{1}){i,2} + pi;
    else means.(m{1}).(q{1}){i,2} =  means.(m{1}).(q{1}){i,2};
     end
    
      if means.(m{1}).(q{1}){i,3} < 0
        means.(m{1}).(q{1}){i,3} =  means.(m{1}).(q{1}){i,3} + pi;
    else means.(m{1}).(q{1}){i,3} =  means.(m{1}).(q{1}){i,3};
      end
    
       if means.(m{1}).(q{1}){i,4} < 0
        means.(m{1}).(q{1}){i,4} =  means.(m{1}).(q{1}){i,4} + pi;
    else means.(m{1}).(q{1}){i,4} =  means.(m{1}).(q{1}){i,4};
       end
    
       
       if means.(m{1}).(q{1}){i,5} < 0
        means.(m{1}).(q{1}){i,5} =  means.(m{1}).(q{1}){i,5} + pi;
    else means.(m{1}).(q{1}){i,5} =  means.(m{1}).(q{1}){i,5};
    end
        end
  end
end

%jetzt v tests fürt means
%ich loope glaube nicht nach zeit.
%ALL
move = {'First', 'Middle', 'Last'};    
qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
for q = qua
    
        for m = move
vcell = means.(m{1}).(q{1})(:,5) %je nachdem welchen zeitpunkt ich haben will
vvektor = cell2mat(vcell)
[p, v] = circ_vtest([vvektor*2], 90*pi/90) %mal zwei ist die axiale korrektur. in radians ist es eh schon

 vresults.(m{1}).(q{1}).all = {v, p}
end
end

   %FIRSSTQUARTER
   
move = {'First', 'Middle', 'Last'};    
qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
for q = qua
        for m = move
vcell = means.(m{1}).(q{1})(:,1) %je nachdem welchen zeitpunkt ich haben will
vvektor = cell2mat(vcell)
[p, v] = circ_vtest([vvektor*2], 90*pi/90) %mal zwei ist die axiale korrektur. in radians ist es eh schon

 vresults.(m{1}).(q{1}).erstesviertel = {v, p}
end
end

      %SECONDQUARTER

move = {'First', 'Middle', 'Last'};    
qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
for q = qua
        for m = move
vcell = means.(m{1}).(q{1})(:,2) %je nachdem welchen zeitpunkt ich haben will
vvektor = cell2mat(vcell)
[p, v] = circ_vtest([vvektor*2], 90*pi/90) %mal zwei ist die axiale korrektur. in radians ist es eh schon

 vresults.(m{1}).(q{1}).zweitesviertel = {v, p}
end
end
        
%THIRDQUARTER
      
move = {'First', 'Middle', 'Last'};    
qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
for q = qua
        for m = move
vcell = means.(m{1}).(q{1})(:,3) %je nachdem welchen zeitpunkt ich haben will
vvektor = cell2mat(vcell)
[p, v] = circ_vtest([vvektor*2], 90*pi/90) %mal zwei ist die axiale korrektur. in radians ist es eh schon

 vresults.(m{1}).(q{1}).drittesviertel = {v, p}
end
end
        
      %FOURTHQUARTER
      
move = {'First', 'Middle', 'Last'};    
qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
for q = qua
        for m = move
vcell = means.(m{1}).(q{1})(:,4) %je nachdem welchen zeitpunkt ich haben will
vvektor = cell2mat(vcell)
[p, v] = circ_vtest([vvektor*2], 90*pi/90) %mal zwei ist die axiale korrektur. in radians ist es eh schon

 vresults.(m{1}).(q{1}).viertessviertel = {v, p}
end
end


%jetzt histogrammme
%hierfür sollte es von radians wieder in angle transformiert werden
      %das funzt so. nein ERSST QUALI MAANANANNAN
%move = {'First', 'Middle', 'Last'};    
%qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
%for q = qua
  %      for m = move
%datahisto = cell2mat(means.(m{1}).(q{1}))
%histotesti = circ_rad2ang(datahisto(:,5))
%histogrami = CircHist(histotesti, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true)

%toPng([histogrami])

%close all
%close all hidden  
%        end
%end

%interessanter: histogram nochmal mit allen werten... aussagekräöftiger
%warhsch...

%move = {'First', 'Middle', 'Last'};    
%qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
%for q = qua
 %       for m = move
%vektorall = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})']


%histogrami = CircHist(vektorall, [0:10:180], 'areAxialData', true, 'pointReflectAxialData', true)


%toPng([histogrami])

%close all
%close all hidden  
 %       end
%end
 
%mean stroke anzahl for single subjects
 for vp = 1:16
     meanstroke.subj(vp) = mean(expdata.subj(vp).numstrokes)
 end
 
 %anova STROKE

numstrokesall = [expdata.subj(1).numstrokes; expdata.subj(2).numstrokes; expdata.subj(3).numstrokes; expdata.subj(4).numstrokes; expdata.subj(5).numstrokes; expdata.subj(6).numstrokes; expdata.subj(7).numstrokes; expdata.subj(8).numstrokes; expdata.subj(9).numstrokes; expdata.subj(10).numstrokes; expdata.subj(11).numstrokes; expdata.subj(12).numstrokes; expdata.subj(13).numstrokes; expdata.subj(14).numstrokes; expdata.subj(15).numstrokes; expdata.subj(16).numstrokes]
qualityall = [Results_VP1(:,2); Results_VP2(:,2); Results_VP3(:,2); Results_VP4(:,2); Results_VP5(:,2); Results_VP6(:,2); Results_VP7(:,2); Results_VP8(:,2); Results_VP9(:,2); Results_VP10(:,2); Results_VP11(:,2); Results_VP12(:,2); Results_VP13(:,2); Results_VP14(:,2); Results_VP15(:,2); Results_VP16(:,2);];
tablestrokequa = [numstrokesall, qualityall];
%strokequa = grpstats(tablestrokequa, qualityall,{'mean', 'std'}) %noch sd etc
%[p,tbl,stats] = anova1(numstrokesall, qualityall) %nicht sig %MUSS NOCH REPATED MEASURES WERDEN


%t test frequency unterschied      

for  vp = 1:16
     for trial=1:72     
 if expdata.subj(vp).shifted.zero.First(trial)>=75 && expdata.subj(vp).shifted.zero.First(trial)<=105 
     trialOrthzero(vp,trial)=1;
 else trialOrthzero(vp,trial)=0;
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
 
 
 if expdata.subj(vp).shifted.fifty.First(trial)>=75 && expdata.subj(vp).shifted.fifty.First(trial)<=105 
     trialOrthfifty(vp,trial)=1;
 else trialOrthfifty(vp,trial)=0;
 end
  
     end
   
end

for i = 1:16
prozentvektorfirstzero(i, 1) = mean(trialOrthzero(i,:));
prozentvektorfirstfifty(i, 1) = mean(trialOrthfifty(i,:));
end

% T TEST
mean(prozentvektorfirstzero)
mean(prozentvektorfirstfifty)
[t,p,ci,stats] = ttest(prozentvektorfirstzero, prozentvektorfirstfifty)


%v test nicht means sondern lle daten


move = {'First', 'Middle', 'Last'};    
qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
for q = qua
        for m = move
vektorall1 = [expdata.subj(1).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(2).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(3).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(4).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(5).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(6).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(7).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(8).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(9).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(10).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(11).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(12).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(13).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(14).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(15).shifted.(q{1}).(m{1})(1,1:17)'; expdata.subj(16).shifted.(q{1}).(m{1})(1,1:17)'];
vektorall2 = [expdata.subj(1).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(2).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(3).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(4).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(5).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(6).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(7).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(8).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(9).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(10).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(11).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(12).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(13).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(14).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(15).shifted.(q{1}).(m{1})(1,18:36)'; expdata.subj(16).shifted.(q{1}).(m{1})(1,18:36)'];
vektorall3 = [expdata.subj(1).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(2).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(3).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(4).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(5).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(6).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(7).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(8).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(9).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(10).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(11).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(12).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(13).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(14).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(15).shifted.(q{1}).(m{1})(1,37:55)'; expdata.subj(16).shifted.(q{1}).(m{1})(1,37:55)'];
vektorall4 = [expdata.subj(1).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(2).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(3).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(4).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(5).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(6).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(7).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(8).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(9).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(10).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(11).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(12).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(13).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(14).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(15).shifted.(q{1}).(m{1})(1,56:72)'; expdata.subj(16).shifted.(q{1}).(m{1})(1,56:72)'];
vektorallall = [expdata.subj(1).shifted.(q{1}).(m{1})'; expdata.subj(2).shifted.(q{1}).(m{1})'; expdata.subj(3).shifted.(q{1}).(m{1})'; expdata.subj(4).shifted.(q{1}).(m{1})'; expdata.subj(5).shifted.(q{1}).(m{1})'; expdata.subj(6).shifted.(q{1}).(m{1})'; expdata.subj(7).shifted.(q{1}).(m{1})'; expdata.subj(8).shifted.(q{1}).(m{1})'; expdata.subj(9).shifted.(q{1}).(m{1})'; expdata.subj(10).shifted.(q{1}).(m{1})'; expdata.subj(11).shifted.(q{1}).(m{1})'; expdata.subj(12).shifted.(q{1}).(m{1})'; expdata.subj(13).shifted.(q{1}).(m{1})'; expdata.subj(14).shifted.(q{1}).(m{1})'; expdata.subj(15).shifted.(q{1}).(m{1})'; expdata.subj(16).shifted.(q{1}).(m{1})'];

vektorall1 = circ_ang2rad(vektorall1);
vektorall2 = circ_ang2rad(vektorall2);
vektorall3 = circ_ang2rad(vektorall3);
vektorall4 = circ_ang2rad(vektorall4);
vektorallall = circ_ang2rad(vektorallall); 

[p, v] = circ_vtest([vektorall1*2], 90*pi/90)
  vresultsall.(m{1}).(q{1}).erstesviertel = {v, p};
  
  
[p, v] = circ_vtest([vektorall2*2], 90*pi/90)
  vresultsall.(m{1}).(q{1}).zweitesviertel = {v, p};
  
  
[p, v] = circ_vtest([vektorall3*2], 90*pi/90)
  vresultsall.(m{1}).(q{1}).drittessviertel = {v, p};
  
  
[p, v] = circ_vtest([vektorall4*2], 90*pi/90)
  vresultsall.(m{1}).(q{1}).viertesviertel = {v, p};
  
  
[p, v] = circ_vtest([vektorallall*2], 90*pi/90)
  vresultsall.(m{1}).(q{1}).allall = {v, p};
        end
end



%descriptives



mean(numstrokesall)
std(numstrokesall)

timeall = [expdata.subj(1).timeonstim; expdata.subj(2).timeonstim; expdata.subj(3).timeonstim; expdata.subj(4).timeonstim; expdata.subj(5).timeonstim; expdata.subj(6).timeonstim; expdata.subj(7).timeonstim; expdata.subj(8).timeonstim; expdata.subj(9).timeonstim; expdata.subj(10).timeonstim; expdata.subj(11).timeonstim; expdata.subj(12).timeonstim; expdata.subj(13).timeonstim; expdata.subj(14).timeonstim; expdata.subj(15).timeonstim; expdata.subj(16).timeonstim]

mean(timeall)
std(timeall)


changesall = [0]
for i = 1:16
    changesall = [changesall; expdata.subj(i).nrchanges'];
end
changesall(1) = [];
mean(changesall)
std(changesall)


accuracyall = [0]
for i = 1:16
   accuracyall = [accuracyall; expdata.subj(i).accuracy];
end
accuracyall(1) = [];
mean(accuracyall)
std(accuracyall)

%VLLT AUCH NOCH V TEST FÜR JEDE PERSON EINZELN ja ist in "weiteres ganz
%unten