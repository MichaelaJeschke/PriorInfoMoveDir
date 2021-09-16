%für first mmiddle last ZERO wegen replikation nochmal przentortho und
%rlength anova. dafür je 3 spalten erstellen für spss mit vp pro zeile
%means...

%orthovektr erstellen
for  vp = 1:16
     for trial=1:72     
 if expdata.subj(vp).shifted.zero.First(trial)>=75 && expdata.subj(vp).shifted.zero.First(trial)<=105 
     trialOrthzerofirst(vp,trial)=1;
 else trialOrthzerofirst(vp,trial)=0;
 end
     end
end
for  vp = 1:16
     for trial=1:72     
 if expdata.subj(vp).shifted.zero.Middle(trial)>=75 && expdata.subj(vp).shifted.zero.Middle(trial)<=105 
     trialOrthzeromiddle(vp,trial)=1;
 else trialOrthzeromiddle(vp,trial)=0;
 end
     end
end
for  vp = 1:16
     for trial=1:72     
 if expdata.subj(vp).shifted.zero.Last(trial)>=75 && expdata.subj(vp).shifted.zero.Last(trial)<=105 
     trialOrthzerolast(vp,trial)=1;
 else trialOrthzerolast(vp,trial)=0;
 end
     end
end

for i = 1:16
prozentvektorzerofirst(i, 1) = mean(trialOrthzerofirst(i,:));
prozentvektorzeromiddle(i, 1) = mean(trialOrthzeromiddle(i,:));
prozentvektorzerolast(i, 1) = mean(trialOrthzerolast(i,:));
end

anovafirstmiddlelastmatrix = [prozentvektorzerofirst, prozentvektorzeromiddle, prozentvektorzerolast], 





%varianzüberprüfung erstmal array mit allen r-werten für qualities first
%move für jede vp
%da ich das mit spss machen werde mach 5 spalten für jede quaity mit den r
%werten für die vps. für first move nur
 move = {'First', 'Middle', 'Last'};    
qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
for i= 1:16
for q = qua
        for m = move
    radiansvek = circ_ang2rad(expdata.subj(i).shifted.(q{1}).(m{1}))';
rlength.(m{1}).(q{1}){i} = circ_r(radiansvek*2);
        end
end
end
%MUSS GUCKEN OB DIE MAL ZWEI HIER AUCH DIE RICHTIGE KORREKTUR FRÜRAXIAL
%IST::::::

%NUN ANOVA VGL ZWISCHENDEN QUALITIES matrix erstelen
zeroR = cell2mat(rlength.First.zero)';
fifteenR= cell2mat(rlength.First.fifteen)';
twentyfiveR= cell2mat(rlength.First.twentyfive)';
thirtyfiveR= cell2mat(rlength.First.thirtyfive)';
fiftyR = cell2mat(rlength.First.fifty)';

anovaspssR = [zeroR fifteenR twentyfiveR thirtyfiveR fiftyR];

[t,p,ci,stats] = ttest(zeroR, fiftyR) %der wird sig die anova nicht


%gucken: gibt es profesionallen test für varianz? gibt es einen profitest für
%paired mittelwert? 


%rlength first middle last ANOVA

vektorrlengthfirst = cell2mat(rlength.First.zero);
vektorrlengthmiddle = cell2mat(rlength.Middle.zero);
vektorrlengthlast = cell2mat(rlength.Last.zero);
anovafirstmiddlelastrlength = [vektorrlengthfirst', vektorrlengthmiddle', vektorrlengthlast']


%quality


p = [1:16];
vp = 0;
for person = p
    vp = vp+1;
    
    for i = 1:360
        %FIRST %quality 0 cell
        if expdata.subj(vp).quality(i) == 0
        expdata.subj(vp).accuracies.zero(i) = expdata.subj(vp).accuracy(i);
        else expdata.subj(vp).accuracies.zero(i)  = NaN;
        end   
        %15   
        if expdata.subj(vp).quality(i) == 15
        expdata.subj(vp).accuracies.fifteen(i) =  expdata.subj(vp).accuracy(i);
        else expdata.subj(vp).accuracies.fifteen(i)  = NaN;
        end
        %25
        if expdata.subj(vp).quality(i) == 25
        expdata.subj(vp).accuracies.twentyfive(i) =  expdata.subj(vp).accuracy(i);
        else expdata.subj(vp).accuracies.twentyfive(i)  = NaN;
        end
         %35
         if expdata.subj(vp).quality(i) == 35
         expdata.subj(vp).accuracies.thirtyfive(i) =  expdata.subj(vp).accuracy(i);
         else expdata.subj(vp).accuracies.thirtyfive(i)  = NaN;
         end
         %50
         if expdata.subj(vp).quality(i) == 50
         expdata.subj(vp).accuracies.fifty(i) =  expdata.subj(vp).accuracy(i);
         else expdata.subj(vp).accuracies.fifty(i)  = NaN;
         end
    end

        
expdata.subj(vp).accuracies.zero(isnan(expdata.subj(vp).accuracies.zero)) = [];
expdata.subj(vp).accuracies.fifteen(isnan( expdata.subj(vp).accuracies.fifteen)) = [];
expdata.subj(vp).accuracies.twentyfive(isnan(expdata.subj(vp).accuracies.twentyfive)) = [];
expdata.subj(vp).accuracies.thirtyfive(isnan(expdata.subj(vp).accuracies.thirtyfive)) = [];
expdata.subj(vp).accuracies.fifty(isnan(expdata.subj(vp).accuracies.fifty)) = [];

end

%KATASTROPHE
% ERST MEAN VEKTOR MACHEN MIT ALLEN 16 MEAN ACCURAIES für jede quality und dann auf jeden
% wert die transformation anwenden

%MUSS NOCH VEKTOREN EINZELNE FÜR QAULITES MIT DEREN ACCURACY WERTE MACHEN
%DIE ICH DANN EINSETZTEN KANN FÜR MEAN BERECHNUNG
for i = 1:16
    
    meanaccuracy(i,1) = mean(expdata.subj(i).accuracies.zero);
    meanaccuracy(i,2) = mean(expdata.subj(i).accuracies.fifteen);
     meanaccuracy(i,3) = mean(expdata.subj(i).accuracies.twentyfive);
      meanaccuracy(i,4) = mean(expdata.subj(i).accuracies.thirtyfive);
       meanaccuracy(i,5) = mean(expdata.subj(i).accuracies.fifty);
end

meanaccuracy_rau = rau(meanaccuracy);
%***********************BIS HIER OK ERSTMAL MATRIX KURZ SPEICHERN UND MIT
%SPSS CHECKEN

writematrix(meanaccuracy_rau,'accuracyagain.csv') 

for i = 1:5760
    if qualityall(i) == 0
        quality2(i) = 1;
    end
    if qualityall(i) == 15
        quality2(i) = 2;
    end
     if qualityall(i) == 25
        quality2(i) = 3;
     end
      if qualityall(i) == 35
        quality2(i) = 4;
      end
      if qualityall(i) == 50
        quality2(i) = 5;
    end
        
end

quality2 = quality2'; 

subjid3= 0;
for i = 1:16
subject = repmat((i), 360, 1);
subjid3 = [subjid3; subject];
end
subjid3(1) = [];

%anovamatrixacc = [accuracyrau  quality2 subjid3];
%[RMAOV1] = RMAOV1(anovamatrixacc)
%ff 15  4    0.000 p 1.0000
%kommt was komisches raus trotz rau transformation ja weils FLAASCH WAR....

%texture
%ersmtal stimpair für alle brauchbar machen
stimpairall = [0 0]
for i = 1:16
stimpairall = [stimpairall; expdata.subj(i).stimpair];

end
stimpairall(1, :) = [];

for i = 1:5760
    if stimpairall((i), 1) == 21 && stimpairall((i), 2) == 31
        stimpairallnew(i) = 1;
    end
    if stimpairall((i), 1) == 31 && stimpairall((i), 2) == 21
        stimpairallnew(i) = 1;
    end
    
    if stimpairall((i), 1) == 31 && stimpairall((i), 2) == 41
        stimpairallnew(i) = 2;
    end
    if stimpairall((i), 1) == 41 && stimpairall((i), 2) == 31
        stimpairallnew(i) = 2;
    end
    
    if stimpairall((i), 1) == 51 && stimpairall((i), 2) == 41
        stimpairallnew(i) = 3;
    end
    if stimpairall((i), 1) == 41 && stimpairall((i), 2) == 51
        stimpairallnew(i) = 3;
    end
end
stimpairallnew = stimpairallnew';

%anovamatrixtex = [accuracyall stimpairallnew subjid3];

%[RMAOV1] = RMAOV1(anovamatrixtex)
%F 15, 2  0.101   p .904

%vllt mal rechnen ob überzufällig viele in orthogonalem fenster liegen bei
%zero First 
%t test against 16.6



%v test für einzelne

 move = {'First', 'Middle', 'Last'};    
qua = {'zero','fifteen','twentyfive', 'thirtyfive', 'fifty'}; 
for i= 1:16
for q = qua
        for m = move
            
vcell = expdata.subj(i).shifted.(q{1}).(m{1});
vvektor = circ_ang2rad(vcell);
[p, v] = circ_vtest([vvektor*2], 90*pi/90); %bei 90 rechne ich nbur durch 90 bei umwaldnlung in rad wegen axialer korrektur

 vresultssingle.(m{1}).(q{1}){i} = {v,p};
 
        end
end
end


%anova zeit


timeall = [expdata.subj(1).timeonstim; expdata.subj(2).timeonstim; expdata.subj(3).timeonstim; expdata.subj(4).timeonstim; expdata.subj(5).timeonstim; expdata.subj(6).timeonstim; expdata.subj(7).timeonstim; expdata.subj(8).timeonstim; expdata.subj(9).timeonstim; expdata.subj(10).timeonstim; expdata.subj(11).timeonstim; expdata.subj(12).timeonstim; expdata.subj(13).timeonstim; expdata.subj(14).timeonstim; expdata.subj(15).timeonstim; expdata.subj(16).timeonstim]
qualityall = [Results_VP1(:,2); Results_VP2(:,2); Results_VP3(:,2); Results_VP4(:,2); Results_VP5(:,2); Results_VP6(:,2); Results_VP7(:,2); Results_VP8(:,2); Results_VP9(:,2); Results_VP10(:,2); Results_VP11(:,2); Results_VP12(:,2); Results_VP13(:,2); Results_VP14(:,2); Results_VP15(:,2); Results_VP16(:,2);];
%tabletimequa = [timeall, quality2, subjid3];
%timequa = grpstats(tabletimequa, quality2,{'mean'}) %noch sd etc
%[p,tbl,stats] = anova1(timeall, qualityall) %nicht sig %MUSS NOCH REPATED MEASURES WERDEN
%[RMAOV1] = RMAOV1(tabletimequa)

%t test agianst 16.6 für zero orthis
[h,p] = ttest(prozentvektorfirstzero,16.67)



