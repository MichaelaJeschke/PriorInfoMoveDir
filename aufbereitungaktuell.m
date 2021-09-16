p = [1:16];
vp = 0;
for person = p
    vp = vp+1;
load(['Results/S' num2str(person) '_Exploration']);
load(['Results/S' num2str(person) '_Results']);
load(['Results/S' num2str(person) '_Stim']);
load(['Results/S' num2str(person) '_timeANDanswer']);
load(['Results/S' num2str(person) '_Moves']);


expdata.subj(vp).nrchanges = zeros(1,360);
expdata.subj(vp).nrchanges = no_cha;
expdata.subj(vp).angles.comp = Angle.Comp;
expdata.subj(vp).angles.ref = Angle.Ref;
expdata.subj(vp).angles.tex = Angle.texture;
expdata.subj(vp).moves.First = move.first;
expdata.subj(vp).moves.Middle = move.middle;
expdata.subj(vp).moves.Last = move.last;

%number of strokes added
for i = 1:360
expdata.subj(vp).numstrokes(i,1) = length(expdata.subj(vp).angles.ref {1,i}) + length(expdata.subj(vp).angles.comp{1,i});
end

%time on stim1 (der erstgewählte) und tinme on stim2 addiert
for i = 1:360
expdata.subj(vp).timeonstim(i,1) = time_on_stim1(1,i) +  time_on_stim2(1,i);
end

%accuracy
for i = 1:360
expdata.subj(vp).accuracy(i,1) = responses.cor(1,i);
end

%stimpair evtl für frequency effect und um erst und zweigewählter
%auseinander zu halten oder so. first column is first touched. smaller value one is comp
expdata.subj(vp).stimpair = StimPair;



%now: moves transformieren zu 180-0 statt 0   -90/90   0

for i = 1:360
    if expdata.subj(vp).moves.First(i) < 0
        expdata.subj(vp).moves.First(i) =  expdata.subj(vp).moves.First(i) + 180;
    else expdata.subj(vp).moves.First(i) =  expdata.subj(vp).moves.First(i);
    end
end


for i = 1:360
    if expdata.subj(vp).moves.Middle(i) < 0
        expdata.subj(vp).moves.Middle(i) =  expdata.subj(vp).moves.Middle(i) + 180;
    else expdata.subj(vp).moves.Middle(i) =  expdata.subj(vp).moves.Middle(i);
    end
end


for i = 1:360
    if expdata.subj(vp).moves.Last(i) < 0
        expdata.subj(vp).moves.Last(i) =  expdata.subj(vp).moves.Last(i) + 180;
    else expdata.subj(vp).moves.Last(i) =  expdata.subj(vp).moves.Last(i);
    end
end


%now: realignement für the histograms. die texture der x achse anspassen

%ACHTUNG NOCH DIE EINSEN FILTER BEFORE SHIFTING die gehen sosnt verloren!!!

%und die moves demenstrechend shiften.
%FIRST
for i = 1:360
    orientation = cell2mat(expdata.subj(vp).angles.tex(i));
shiftID = 360 - orientation;
expdata.subj(vp).shifted.First(i) = expdata.subj(vp).moves.First(i) + (shiftID); 


%für 180:360 minus 180
%für > 360 minus 360

if expdata.subj(vp).shifted.First(i) > 180 &&  expdata.subj(vp).shifted.First(i) <360
    expdata.subj(vp).shifted.First(i) = expdata.subj(vp).shifted.First(i) - 180;
else expdata.subj(vp).shifted.First(i)=expdata.subj(vp).shifted.First(i);
end

if expdata.subj(vp).shifted.First(i) > 360
    expdata.subj(vp).shifted.First(i) = expdata.subj(vp).shifted.First(i) - 360;
else expdata.subj(vp).shifted.First(i)= expdata.subj(vp).shifted.First(i);
end
end

%MIDDLE
for i = 1:360
    orientation = cell2mat(expdata.subj(vp).angles.tex(i));
shiftID = 360 - orientation;
expdata.subj(vp).shifted.Middle(i) = expdata.subj(vp).moves.Middle(i) + (shiftID); 


%für 180:360 minus 180
%für > 360 minus 360

if expdata.subj(vp).shifted.Middle(i) > 180 &&  expdata.subj(vp).shifted.Middle(i) <360
    expdata.subj(vp).shifted.Middle(i) = expdata.subj(vp).shifted.Middle(i) - 180;
else expdata.subj(vp).shifted.Middle(i)=expdata.subj(vp).shifted.Middle(i);
end

if expdata.subj(vp).shifted.Middle(i) > 360
    expdata.subj(vp).shifted.Middle(i) = expdata.subj(vp).shifted.Middle(i) - 360;
else expdata.subj(vp).shifted.Middle(i)= expdata.subj(vp).shifted.Middle(i);
end

end

%LAST
for i = 1:360
    orientation = cell2mat(expdata.subj(vp).angles.tex(i));
shiftID = 360 - orientation;
expdata.subj(vp).shifted.Last(i) = expdata.subj(vp).moves.Last(i) + (shiftID); 


%für 180:360 minus 180
%für > 360 minus 360

if expdata.subj(vp).shifted.Last(i) > 180 &&  expdata.subj(vp).shifted.Last(i) <360
    expdata.subj(vp).shifted.Last(i) = expdata.subj(vp).shifted.Last(i) - 180;
else expdata.subj(vp).shifted.Last(i)=expdata.subj(vp).shifted.Last(i);
end

if expdata.subj(vp).shifted.Last(i) > 360
    expdata.subj(vp).shifted.Last(i) = expdata.subj(vp).shifted.Last(i) - 360;
else expdata.subj(vp).shifted.Last(i)= expdata.subj(vp).shifted.Last(i);
end

end

end


%quality (leider etwas umständlich aber geht nicht anders keine lust es im
%vorskript zu ändern)


expdata.subj(1).quality = Results_VP1(:,2); 
expdata.subj(2).quality = Results_VP2(:,2); 
expdata.subj(3).quality = Results_VP3(:,2); 
expdata.subj(4).quality = Results_VP4(:,2); 
expdata.subj(5).quality = Results_VP5(:,2); 
expdata.subj(6).quality = Results_VP6(:,2); 
expdata.subj(7).quality = Results_VP7(:,2); 
expdata.subj(8).quality = Results_VP8(:,2); 
expdata.subj(9).quality = Results_VP9(:,2); 
expdata.subj(10).quality = Results_VP10(:,2); 
expdata.subj(11).quality = Results_VP11(:,2);
expdata.subj(12).quality = Results_VP12(:,2);
expdata.subj(13).quality = Results_VP13(:,2);
expdata.subj(14).quality = Results_VP14(:,2);
expdata.subj(15).quality = Results_VP15(:,2);
expdata.subj(16).quality = Results_VP16(:,2);



%RAUSFINDEN WIE MAN NUN NACH QUALITY AUFTEILT:::FÜR HISTOS ABER AUCH SO



p = [1:16];
vp = 0;
for person = p
    vp = vp+1;
    
    for i = 1:360
        %FIRST %quality 0 cell
        if expdata.subj(vp).quality(i) == 0
        expdata.subj(vp).shifted.zero.First(i) = expdata.subj(vp).shifted.First(i);
        else expdata.subj(vp).shifted.zero.First(i)  = NaN;
        end   
        %15   
        if expdata.subj(vp).quality(i) == 15
        expdata.subj(vp).shifted.fifteen.First(i) = expdata.subj(vp).shifted.First(i);
        else expdata.subj(vp).shifted.fifteen.First(i)  = NaN;
        end
        %25
        if expdata.subj(vp).quality(i) == 25
        expdata.subj(vp).shifted.twentyfive.First(i) = expdata.subj(vp).shifted.First(i);
        else expdata.subj(vp).shifted.twentyfive.First(i)  = NaN;
        end
         %35
         if expdata.subj(vp).quality(i) == 35
         expdata.subj(vp).shifted.thirtyfive.First(i) = expdata.subj(vp).shifted.First(i);
        else expdata.subj(vp).shifted.thirtyfive.First(i)  = NaN;
         end
         %50
         if expdata.subj(vp).quality(i) == 50
         expdata.subj(vp).shifted.fifty.First(i) = expdata.subj(vp).shifted.First(i);
         else expdata.subj(vp).shifted.fifty.First(i)  = NaN;
         end
         
       %MIDDLE  %quality 0 cell
        if expdata.subj(vp).quality(i) == 0
        expdata.subj(vp).shifted.zero.Middle(i) = expdata.subj(vp).shifted.Middle(i);
        else expdata.subj(vp).shifted.zero.Middle(i)  = NaN;
        end   
        %15   
        if expdata.subj(vp).quality(i) == 15
        expdata.subj(vp).shifted.fifteen.Middle(i) = expdata.subj(vp).shifted.Middle(i);
        else expdata.subj(vp).shifted.fifteen.Middle(i)  = NaN;
        end
        %25
        if expdata.subj(vp).quality(i) == 25
        expdata.subj(vp).shifted.twentyfive.Middle(i) = expdata.subj(vp).shifted.Middle(i);
        else expdata.subj(vp).shifted.twentyfive.Middle(i)  = NaN;
        end
         %35
         if expdata.subj(vp).quality(i) == 35
         expdata.subj(vp).shifted.thirtyfive.Middle(i) = expdata.subj(vp).shifted.Middle(i);
        else expdata.subj(vp).shifted.thirtyfive.Middle(i)  = NaN;
         end
         %50
         if expdata.subj(vp).quality(i) == 50
         expdata.subj(vp).shifted.fifty.Middle(i) = expdata.subj(vp).shifted.Middle(i);
         else expdata.subj(vp).shifted.fifty.Middle(i)  = NaN;
         end
         
        %LAST  %quality 0 cell
        if expdata.subj(vp).quality(i) == 0
        expdata.subj(vp).shifted.zero.Last(i) = expdata.subj(vp).shifted.Last(i);
        else expdata.subj(vp).shifted.zero.Last(i)  = NaN;
        end   
        %15   
        if expdata.subj(vp).quality(i) == 15
        expdata.subj(vp).shifted.fifteen.Last(i) = expdata.subj(vp).shifted.Last(i);
        else expdata.subj(vp).shifted.fifteen.Last(i)  = NaN;
        end
        %25
        if expdata.subj(vp).quality(i) == 25
        expdata.subj(vp).shifted.twentyfive.Last(i) = expdata.subj(vp).shifted.Last(i);
        else expdata.subj(vp).shifted.twentyfive.Last(i)  = NaN;
        end
         %35
         if expdata.subj(vp).quality(i) == 35
         expdata.subj(vp).shifted.thirtyfive.Last(i) = expdata.subj(vp).shifted.Last(i);
        else expdata.subj(vp).shifted.thirtyfive.Last(i)  = NaN;
         end
         %50
         if expdata.subj(vp).quality(i) == 50
         expdata.subj(vp).shifted.fifty.Last(i) = expdata.subj(vp).shifted.Last(i);
         else expdata.subj(vp).shifted.fifty.Last(i)  = NaN;
         end
    end
end

%NOCH ALLE NAN RAUSWERFEN
p = [1:16];
vp = 0;
for person = p
    vp = vp+1;
    
expdata.subj(vp).shifted.zero.First(isnan(expdata.subj(vp).shifted.zero.First)) = [];
expdata.subj(vp).shifted.fifteen.First(isnan(expdata.subj(vp).shifted.fifteen.First)) = [];
expdata.subj(vp).shifted.twentyfive.First(isnan(expdata.subj(vp).shifted.twentyfive.First)) = [];
expdata.subj(vp).shifted.thirtyfive.First(isnan(expdata.subj(vp).shifted.thirtyfive.First)) = [];
expdata.subj(vp).shifted.fifty.First(isnan(expdata.subj(vp).shifted.fifty.First)) = [];
expdata.subj(vp).shifted.zero.Middle(isnan(expdata.subj(vp).shifted.zero.Middle)) = [];
expdata.subj(vp).shifted.fifteen.Middle(isnan(expdata.subj(vp).shifted.fifteen.Middle)) = [];
expdata.subj(vp).shifted.twentyfive.Middle(isnan(expdata.subj(vp).shifted.twentyfive.Middle)) = [];
expdata.subj(vp).shifted.thirtyfive.Middle(isnan(expdata.subj(vp).shifted.thirtyfive.Middle)) = [];
expdata.subj(vp).shifted.fifty.Middle(isnan(expdata.subj(vp).shifted.fifty.Middle)) = [];
expdata.subj(vp).shifted.zero.Last(isnan(expdata.subj(vp).shifted.zero.Last)) = [];
expdata.subj(vp).shifted.fifteen.Last(isnan(expdata.subj(vp).shifted.fifteen.Last)) = [];
expdata.subj(vp).shifted.twentyfive.Last(isnan(expdata.subj(vp).shifted.twentyfive.Last)) = [];
expdata.subj(vp).shifted.thirtyfive.Last(isnan(expdata.subj(vp).shifted.thirtyfive.Last)) = [];
expdata.subj(vp).shifted.fifty.Last(isnan(expdata.subj(vp).shifted.fifty.Last)) = [];
end


%run rumgeteste




