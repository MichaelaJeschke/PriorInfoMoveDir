% Auswertung pro VP 
clear all
close all

Subjects = 16; % Number of participants
%sourcePath = 'C:\Users\Aaron\Desktop\Paper_und_Projekt\Projekt 5 PriorInfoMovementDirection\Experimental Data\Data TestVP\';
sourcePath = 'C:\Users\michi\OneDrive\Desktop\datamichi\'; %hier wird alles hergenommen
savePath = 'C:/Users/michi/OneDrive/Desktop/Results/'; %hierhin werden die sachen gespeichert
trials=1:180; % 360/2 weil 2 blöcke (180) Trials in original, 60 in test
blockind = 1:2; % Number of blocks in original 2 in test 1!
%responses = struct(zeros(length(trials), Subjects)); %% ACHTUNG hier haben
%wir was rausgenommen weil es sonst fehlermeldung gab


for vp=[1:16]; % Number of participants (16 insgesamt, 1:16, hier nur grad 1 zum test)

    
%% Loading the Data (1 at the end is the exp, 0 is the test. sometimes there is another 1 at the beginning this is due to problems with the phantom
    
if vp == 1
    subInitials = 'JN1011';
elseif vp == 2
    subInitials = 'AG1021';
elseif vp == 3
    subInitials = 'JL031';
elseif vp == 4
    subInitials = 'FG1041';
elseif vp == 5
    subInitials = 'LR051';
elseif vp == 6
    subInitials = 'NP061';
elseif vp == 7
    subInitials = 'ST071';
elseif vp == 8
    subInitials = 'SO1081';
elseif vp == 9
    subInitials = 'LS1091';
elseif vp == 10
    subInitials = 'EG101';
elseif vp == 11
    subInitials = 'EM111';
elseif vp == 12
    subInitials = 'MY121';
elseif vp == 13
    subInitials = 'AW1131';
elseif vp == 14
    subInitials = 'LH1131'; %NOCHMAL ÜBERPRÜFEN
elseif vp == 15
    subInitials = 'AM1151';
elseif vp == 16
    subInitials = 'AW161';
end

sourceFilePath = sourcePath; %ist der glieche wie der wos herkommt
sourceFilePath = strcat (sourceFilePath, subInitials); % verkettet strings. bei literal strings anführungszeichen verwenden, hier ohne)
fileName = strcat (sourceFilePath, '.raw'); %siehst du


[   file_blockNr, ...           % BlockNumber
    file_trialNr, ...           % TrialNumber
    file_dataPriorIndex,...     % Image Index given by Programm to select the correct image   
    file_TrialOrientation,...   % [15:30:165]
    file_TrialQuality,...       % Prior Inofrmation Quality
    file_startStim, ...          % Startingposition of stimulus 1eft = (0), rigth = (1)      
    file_refStim, ...           % ref stimulus size
    file_compStim, ...          % comp stimulus size
    file_leftStim, ...          % Stimulus on left position
    file_chosenStim, ...        % chosen Stimulus (die antwort) (daraus wird accuracy gemacht)
    file_weightStimpair ...     % Stimulus weight
] = textread (fileName, ...
    '%f %f %s %f %f %f %f %f %f %f %f');
    fileFull = load(fileName); % (1) BlockNr; (2) TrialNr; (3) PriorIndex; (4) Angle; (5) PriorQuality; (6) LeftStim (left = 0, right = 1); (7) ref stim; (8) comp stim; (9) Left stim (10) Chosen Stim (11) WeightStim
    

tCounter=1; % trial counter
analysedTrialsNr = 1;

for Block = blockind
for trInd=trials
    LEFT  = 0;
    RIGHT = 1;
    REF  = 0;
    COMP = 1;

    ExplError=zeros(1,length(trials));
    
       %% correct answer:

         if (file_refStim(trInd)< file_compStim(trInd)) % bei frq heisst kleinere Zahl höhere Frq. man soll die höhere wählen
            if (file_chosenStim(trInd) == file_refStim(trInd)) %wenn also ref gewählt wurde (kleinere zahl), dann lag man korrekt
                 cor_answ=1; %1 ist korrekt
            else cor_answ=0; %0 ist inkorrekt
            end
         elseif file_compStim(trInd) < file_refStim (trInd) %wenn aber comp die kleinere zahl hat
           if file_chosenStim(trInd) == file_compStim(trInd)%dann sollte auch comp ausgewählt worden sein
                cor_answ=1;
           else cor_answ=0;
           end
         end

         stdChoosen=(file_chosenStim(trInd)==file_refStim(trInd));
         responses.cor(tCounter)=cor_answ; %die hier zeigt für jedes trial ob richtig oder falsch
         responses.std(tCounter)=stdChoosen;
  
%%   StimPair: [1st 2nd] 

% Which one is the left and which one is the right stim
   leftStim= file_leftStim(tCounter);
       if (leftStim == file_refStim(tCounter))
           rightStim = file_compStim(tCounter);
       elseif (leftStim == file_compStim(tCounter))
           rightStim = file_refStim(tCounter);
       end 

% 1eft = (0), rigth = (1) -> starting position (first second stim) 
   if (file_startStim(tCounter) == 1) % Start right, second stim = left
       secondStim = file_leftStim(tCounter);
       if (secondStim == file_refStim(tCounter))
           firstStim = file_compStim(tCounter);
       elseif (secondStim == file_compStim(tCounter))
           firstStim = file_refStim(tCounter);
       end
   elseif (file_startStim(tCounter) == 0)
       firstStim = file_leftStim(tCounter); % start = left
       if (firstStim == file_refStim(tCounter))
           secondStim = file_compStim(tCounter);
       elseif (firstStim == file_compStim(tCounter))
           secondStim = file_refStim(tCounter);
       end
   end
   StimPair(trInd+(180*(Block-1)),:)=[firstStim secondStim]; %hier wird erstgewählter und zweitgewählter aufgelistet

   %%  Einzelne Trials einlesen:

        block = Block-1; % C++ startet bei 0, Matlab bei 1!
        trial = trInd -1; % C++ startet bei 0, Matlab bei 1!
    trjFileName = ([sourceFilePath, '\', 'ET#_', num2str(block), '_', num2str(trial), '_.trj']); %diese ET sind die TRJ dateien

        % Get data from the current .trj file
        [   file_time, ...
            file_xPos, ...
            file_yPos, ...
            file_zPos, ...
            file_Force, ...
            file_aveYVel, ...
            file_stimPos, ... % NO STIM:-1; L0 R1
            file_INstimPos, ... % NOT INstim: -1; OUTERSTIM=1; INNERSTIM=0
        ] = textread (trjFileName, '%f %f %f %f %f %f %f %f');
       
    %% time on stimulus:
index=[];
index=find(file_stimPos~=-1);
        if (file_stimPos(index(1))==0) %1st vs 2nd Stim ALSO ERST UND ZWEITGEWÄHLTER
        time_on_stim1(tCounter)=length(find(file_stimPos==0))*3/(1000);% da messung alle 3 msec geschah --> zeit in sec
        time_on_stim2(tCounter)=length(find(file_stimPos==1))*3/(1000);
        else
        time_on_stim1(tCounter)=length(find(file_stimPos==1))*3/(1000);
        time_on_stim2(tCounter)=length(find(file_stimPos==0))*3/(1000);
        end
        
        if (file_refStim(tCounter)==rightStim) %Ref vs Comp
        time_on_Comp(tCounter)=length(find(file_stimPos==0))*3/(1000);% da messung alle 3 msec geschah --> zeit in sec
        time_on_Ref(tCounter)=length(find(file_stimPos==1))*3/(1000);
        else
        time_on_Comp(tCounter)=length(find(file_stimPos==1))*3/(1000);
        time_on_Ref(tCounter)=length(find(file_stimPos==0))*3/(1000);
        end
     
%% Number of changes
change_ix=[];
change=[];
change_ix=intersect(find(diff(file_stimPos)~=0), find(file_stimPos(1:end-1)~=(-1))); % finde alle indexe die nicht bei (-1) zwischen stimuli anfange und denen ein veräderung in kodierter position folgt
change=change_ix(find(diff(file_stimPos(change_ix))~=0)); %zählt nicht als wechsel wenn letzter wechsel zum gleichen ziel

no_cha(tCounter)=length(change); % Anzahl der Wechsel
       
    
 %% Strokes    --> 1) die von Rand zu Rand

% Eckpunkte (Wechsel von OUTERstim zu INNERstim)
indx_StimAreaStart_uf=sort([intersect(find([diff(file_INstimPos)]==-1),find(file_INstimPos==1)); intersect(find([diff(file_INstimPos)]==1),find(file_INstimPos==-1))]); % Wechsel von Au?enbereich nach Innen
indx_StimAreaEnd_uf=sort([intersect(find([diff(file_INstimPos)]==1),find(file_INstimPos==0)); intersect(find([diff(file_INstimPos)]==-1),find(file_INstimPos==0))]); 
if indx_StimAreaStart_uf(1)>indx_StimAreaEnd_uf(1)
    indx_StimAreaEnd_uf=indx_StimAreaEnd_uf(2:end);
end

indx_StimAreaStart=indx_StimAreaStart_uf(find(indx_StimAreaStart_uf+1~=indx_StimAreaEnd_uf)); % rausfiltern von wacklern
indx_StimAreaEnd=indx_StimAreaEnd_uf(find(indx_StimAreaStart_uf+1~=indx_StimAreaEnd_uf));



%% Strokes (zusaetzliche)--> 2) Richtungsaenderung in X oder Z  

    % Richtungswechsel der Bewegung ohne Beruehrung des Aussenbereichs:
%smoothen:
smooth_x=smoother(file_xPos,20);
smooth_z=smoother(file_zPos,20);

rest=.00000000000001;
slopeX=(diff(smooth_x)./(diff(file_time)+rest));        
slopeZ=(diff(smooth_z)./(diff(file_time)+rest));
slopeZX=(diff(smooth_z)./(diff(smooth_x)));
% knickZX=windowCompareZX(slopeZX, 50);
knicksX=Aaron_windowCompare(slopeX, 50); 
knicksZ=Aaron_windowCompare(slopeZ, 50); 
% 
knicks=union(knicksX,knicksZ); % eins von beiden reicht als Kriterium fuer neuen Stroke
% knicks=knickZX;

%% Force:
    aYForces_row=file_Force-file_weightStimpair(trInd); % Force independend of stimulus-pair weight (Abziehen des GHewichtes des jeweiligen Stimuli-paares)
    aYForces=aYForces_row;
    
        %Ausreisser beseitigen
            % Schritt 1: Vgl mit vorherigen und nach folgenden Wert
                nf=find(aYForces<0);% negative werte raus
                nf=nf(2:(end-1));
                aYForces(nf)=((aYForces(nf-1)+aYForces(nf+1))/2);
   
                
                for l=2:(length(aYForces)-1)
                    if (((aYForces(l-1)+aYForces(l+1))/2) < (aYForces(l)/2))%(aYForces(l)<0)%wenn wert 4mal grösser als drum rum liegende 
                    aYForces(l)=((aYForces(l-1)+aYForces(l+1))/2);
                    end
                end
                
            
            % Schritt 2: Smoothen um rauschen zu erniedrigen
                kernel=15;%
                smoothYForces=smoother(aYForces,kernel);
                touch_indx=find(smoothYForces>=0.1); % mindestent .15N Kraftaufwand

                knicks=intersect(knicks,touch_indx);

% zu nahe Maxima/ doppelte entfernen
    time_intervall=diff(file_time(knicks)); %(2-1)..(3-2)...
     flag=(time_intervall<200);% zwei maxima m?ssen mindestens 200 ms von ein ander entfernt sein
     ind_f=find(flag==1);
     
    
  if (sum(flag)>0 )
     Nknicks=zeros(1,(length(knicks)-sum(flag)));  
     for n=1:sum(flag)
         
        if (smoothYForces(knicks(ind_f(n)))> smoothYForces(knicks(ind_f(n)+1)))
            temp_max_ind(n)=knicks(ind_f(n));
        else temp_max_ind(n)=knicks(ind_f(n)+1);
        end
        
        Nknicks(ind_f(n)-(n-1))=temp_max_ind(n);%
     end
        temp_pair=find(diff(find(flag==1))==1); % sucht den index(innerhalb der zu nahen maximapaare), wenn zwei solche paare hintereinander sind 
         if isempty(temp_pair)==0 % if temp_pair is NOT empty
             
            for tp=1:length(temp_pair)
                if (smoothYForces(temp_max_ind(temp_pair(tp)))> smoothYForces(temp_max_ind(temp_pair(tp)+1))) 
                     Nknicks(ind_f(temp_pair(tp))-(temp_pair(tp)-1))=temp_max_ind(temp_pair(tp));
                else Nknicks(ind_f(temp_pair(tp))-(temp_pair(tp)-1))=temp_max_ind(temp_pair(tp)+1);
                end
            end

         
          end         % else NMaxPressRowIndexes(ind_f(n)-(n-1))=maxPressRowIndexes(ind_f(n)+1);
         zeroPos=find(Nknicks==0);
     
 
     
     z=0;
         for index=1:length(knicks)
             string2eval=[];
             for n=1:sum(flag)
             string2eval=[string2eval '((~(index == ind_f(' num2str(n) '))) & (~(index == ind_f(' num2str(n) ')+1)))' '&'];
             end
             string2eval=string2eval(1:(end-1));
             
             if eval(string2eval)
                 z=z+1;
                Nknicks(zeroPos(z))=knicks(index);
             end 
         end
         
  else
        Nknicks=knicks;
  end
  


    %%

%innerhalb des Ref/Comp Stimulus     
    if (file_refStim(tCounter)==rightStim || file_refStim(tCounter)==rightStim) %% einer der stds rechts
    indx_CompStim= find(file_stimPos==0);
    indx_RefStim= find(file_stimPos==1);
    else indx_CompStim=find(file_stimPos==1);
        indx_RefStim=find(file_stimPos==0);
    end 
    
   %strokes 
    StrokeEnd_Comp=intersect(Nknicks,indx_CompStim);%length(StrokeEnd_Comp)
    StrokeEnd_Ref=intersect(Nknicks,indx_RefStim);%length(StrokeEnd_Ref)
    
    StrokeEnd_Comp=intersect(StrokeEnd_Comp,touch_indx);
    StrokeEnd_Ref=intersect(StrokeEnd_Ref,touch_indx);

%Movement-Start-Point
indx_moveStart_Comp=intersect(indx_StimAreaStart+1,indx_CompStim); %+1 damit nit ein start von nicht-stimulus-bereich zu stimulus-inneren ausfallen w?rde
indx_moveStart_Ref=intersect(indx_StimAreaStart+1,indx_RefStim);
indx_moveEnd_Comp=intersect(indx_StimAreaEnd,indx_CompStim);
indx_moveEnd_Ref=intersect(indx_StimAreaEnd,indx_RefStim);

%% Zus?tzliche Stroke in andere einflechten:
All_indx_moveStart_Comp=indx_moveStart_Comp;
All_indx_moveEnd_Comp=indx_moveEnd_Comp;
All_indx_moveStart_Ref=indx_moveStart_Ref;
All_indx_moveEnd_Ref=indx_moveEnd_Ref;

posPoints_Comp=[];
posPoints_Ref=[];
fin=[];
fin_Ref=[];
ind_short_Ref=[];
ind_short_Ref=[];

if (length(StrokeEnd_Comp)>0)
for i=1:length(StrokeEnd_Comp)
    point=StrokeEnd_Comp(i);
    posPoints_Comp(:,i)=[sum(point>=indx_moveStart_Comp); sum(point>=indx_moveEnd_Comp)];
end
 Within_ind_Comp=(posPoints_Comp(1,:)>posPoints_Comp(2,:))'.*StrokeEnd_Comp;% zwischen einem Start und End
 
 Btw_Comp=[(posPoints_Comp(1,:)==posPoints_Comp(2,:)).* (posPoints_Comp(1,:)~=0)]';% zwischen Ende und n?chstem Start 
 change_ToRight=find(sum([zeros(2,1) diff(posPoints_Comp,[],2)])~=0);
 change_ToLeft=find(fliplr(sum([zeros(2,1) diff(fliplr(posPoints_Comp),[],2)]))~=0);
 posChange_Comp=[intersect(change_ToRight,change_ToLeft) find(change_ToLeft==1)]; % einzelne (nur 1!) punkte zw ende eines stroke und anfang des n?chsten % f?r den allerersten Punkten gilt dass keine ver?nderung von rechts aus ausreicht (da links kein wert)
Btw_Comp(posChange_Comp)=0; % l?schen, wenn nach strokeEnde nur ein Max
 
Btw_ind_Comp=Btw_Comp.*StrokeEnd_Comp;
 
Btw_FirstLast_ind=[intersect(setdiff(change_ToRight,posChange_Comp), find(Btw_ind_Comp)~=0), intersect((ismember(change_ToLeft,change_ToRight)==0).*change_ToLeft, change_ToLeft(2:end))];
Btw_Comp(Btw_FirstLast_ind)=1;
All_indx_moveStart_Comp=sort([All_indx_moveStart_Comp; StrokeEnd_Comp(find(Within_ind_Comp~=0)); StrokeEnd_Comp(intersect(find(Btw_ind_Comp~=1),find(Btw_ind_Comp~=0)))]);
All_indx_moveEnd_Comp=sort([All_indx_moveEnd_Comp; StrokeEnd_Comp(find(Within_ind_Comp~=0)); StrokeEnd_Comp(find(Btw_ind_Comp~=0))]);
end
 
if (length(StrokeEnd_Ref)>0) 
for j=1:length(StrokeEnd_Ref)
    point=StrokeEnd_Ref(j);
    posPoints_Ref(:,j)=[sum(point>=indx_moveStart_Ref); sum(point>=indx_moveEnd_Ref)];
end
 Within_ind_Ref=(posPoints_Ref(1,:)>posPoints_Ref(2,:))'.*StrokeEnd_Ref;

Btw_Ref=[(posPoints_Ref(1,:)==posPoints_Ref(2,:)).* (posPoints_Ref(1,:)~=0)]';% zwischen Ende und n?chstem Start 
 change_ToRight_Ref=find(sum([zeros(2,1) diff(posPoints_Ref,[],2)])~=0);
 change_ToLeft_Ref=find(fliplr(sum([zeros(2,1) diff(fliplr(posPoints_Ref),[],2)]))~=0);
 posChange_Ref=[intersect(change_ToRight_Ref,change_ToLeft_Ref) find(change_ToLeft_Ref==1)]; % einzelne (nur 1!) punkte zw ende eines stroke und anfang des n?chsten % f?r den allerersten Punkten gilt dass keine ver?nderung von rechts aus ausreicht (da links kein wert)
Btw_Ref(posChange_Ref)=0; % l?schen, wenn nach strokeEnde nur ein Max
 
Btw_ind_Ref=Btw_Ref.*StrokeEnd_Ref;
 
Btw_FirstLast_ind_Ref=[intersect((setdiff(change_ToRight_Ref,posChange_Ref)), find(Btw_ind_Ref)~=0), intersect((ismember(change_ToLeft_Ref,change_ToRight_Ref)==0).*change_ToLeft_Ref, change_ToLeft_Ref(2:end))];
Btw_Ref(Btw_FirstLast_ind_Ref)=1;
All_indx_moveStart_Ref=sort([All_indx_moveStart_Ref; StrokeEnd_Ref(find(Within_ind_Ref~=0)); StrokeEnd_Ref(intersect(find(Btw_ind_Ref~=1),find(Btw_ind_Ref~=0)))]);
All_indx_moveEnd_Ref=sort([All_indx_moveEnd_Ref; StrokeEnd_Ref(find(Within_ind_Ref~=0)); StrokeEnd_Ref(find(Btw_ind_Ref~=0))]);
end


%% bereinigen um die Strokes, die zu hoch sind Mittelwertbasiert
Fin_All_indx_moveStart_Ref=[];
Fin_All_indx_moveEnd_Ref=[];
Fin_All_indx_moveStart_Comp=[];
Fin_All_indx_moveEnd_Comp=[];
criterium_yPos=mean(file_yPos(intersect(All_indx_moveStart_Ref(1):All_indx_moveEnd_Ref(end),intersect(find(file_stimPos~=-1),touch_indx))))+2*std(file_yPos(intersect(All_indx_moveStart_Ref(1):All_indx_moveEnd_Ref(end),intersect(find(file_stimPos~=-1),touch_indx))));%+2*std(file_yPos(intersect(All_indx_moveStart_Ref(1):All_indx_moveEnd_Ref(end),intersect(find(file_stimPos~=-1),touch_indx))))%+1/4*std(file_yPos(find(file_stimPos~=-1)))%mean(file_yPos(intersect(find(file_stimPos~=-1),touch_indx)))+1/4*std(file_yPos(find(file_stimPos~=-1)))% f?r vp 2:-10 f?r 1,4,3,10: -8 f?r 10:-6 f?r rest: mean(file_yPos(intersect(find(file_stimPos~=-1),touch_indx)))+1/4*std(file_yPos(find(file_stimPos~=-1)));

% neu in Orientierung
for inl=1:length(All_indx_moveStart_Ref)
    
if ((file_yPos(All_indx_moveStart_Ref(inl))>criterium_yPos) && (file_yPos(All_indx_moveEnd_Ref(inl))<=criterium_yPos))
Fin_All_indx_moveStart_Ref(inl)= [All_indx_moveStart_Ref(inl)+find(file_yPos(All_indx_moveStart_Ref(inl):All_indx_moveEnd_Ref(inl))<criterium_yPos,1)];
Fin_All_indx_moveEnd_Ref(inl)=[All_indx_moveEnd_Ref(inl)];

elseif ((file_yPos(All_indx_moveEnd_Ref(inl))>criterium_yPos) && (file_yPos(All_indx_moveStart_Ref(inl))<=criterium_yPos))
Fin_All_indx_moveStart_Ref(inl)=[All_indx_moveStart_Ref(inl)];  
Fin_All_indx_moveEnd_Ref(inl)=[All_indx_moveStart_Ref(inl)+find(file_yPos(All_indx_moveStart_Ref(inl):All_indx_moveEnd_Ref(inl))>criterium_yPos,1)-1];    

elseif ((file_yPos(All_indx_moveEnd_Ref(inl))>criterium_yPos) && (file_yPos(All_indx_moveStart_Ref(inl))>criterium_yPos))
try
Fin_All_indx_moveStart_Ref(inl)=[All_indx_moveStart_Ref(inl)+find(file_yPos(All_indx_moveStart_Ref(inl):All_indx_moveEnd_Ref(inl))<criterium_yPos,1)];
Fin_All_indx_moveEnd_Ref(inl)=[All_indx_moveStart_Ref(inl)+find(file_yPos(All_indx_moveStart_Ref(inl):All_indx_moveEnd_Ref(inl))>criterium_yPos,1)-1];    
catch
    try
    Fin_All_indx_moveStart_Ref(inl)=[All_indx_moveStart_Ref(inl)+find(file_yPos(All_indx_moveStart_Ref(inl):All_indx_moveEnd_Ref(inl))<criterium_yPos+6,1)];
    Fin_All_indx_moveEnd_Ref(inl)=[All_indx_moveStart_Ref(inl)+find(file_yPos(All_indx_moveStart_Ref(inl):All_indx_moveEnd_Ref(inl))>criterium_yPos+6,1)-1];    
    catch
        Fin_All_indx_moveStart_Ref(inl)=All_indx_moveStart_Ref(inl);
        Fin_All_indx_moveEnd_Ref(inl)=All_indx_moveEnd_Ref(inl);
    end
    
 end
else Fin_All_indx_moveStart_Ref(inl)=All_indx_moveStart_Ref(inl);%(intersect(find(file_yPos(indx_StimAreaStart)<criterium_yPos), find(file_yPos(indx_StimAreaEnd)<criterium_yPos)))%;%;
Fin_All_indx_moveEnd_Ref(inl)=All_indx_moveEnd_Ref(inl);%(intersect(find(file_yPos(indx_StimAreaStart)<criterium_yPos), find(file_yPos(indx_StimAreaEnd)<criterium_yPos)))%All_indx_moveEnd_Ref;%;

end
end


% Comp
criterium_yPos=mean(file_yPos(intersect(All_indx_moveStart_Comp(1):All_indx_moveEnd_Comp(end),intersect(find(file_stimPos~=-1),touch_indx))))+2*std(file_yPos(intersect(All_indx_moveStart_Comp(1):All_indx_moveEnd_Comp(end),intersect(find(file_stimPos~=-1),touch_indx))));%+2*std(file_yPos(intersect(All_indx_moveStart_Ref(1):All_indx_moveEnd_Ref(end),intersect(find(file_stimPos~=-1),touch_indx))))%+1/4*std(file_yPos(find(file_stimPos~=-1)))%mean(file_yPos(intersect(find(file_stimPos~=-1),touch_indx)))+1/4*std(file_yPos(find(file_stimPos~=-1)))% f?r vp 2:-10 f?r 1,4,3,10: -8 f?r 10:-6 f?r rest: mean(file_yPos(intersect(find(file_stimPos~=-1),touch_indx)))+1/4*std(file_yPos(find(file_stimPos~=-1)));

for inl=1:length(All_indx_moveStart_Comp)
    
if ((file_yPos(All_indx_moveStart_Comp(inl))>criterium_yPos) && (file_yPos(All_indx_moveEnd_Comp(inl))<=criterium_yPos))
Fin_All_indx_moveStart_Comp(inl)= [All_indx_moveStart_Comp(inl)+find(file_yPos(All_indx_moveStart_Comp(inl):All_indx_moveEnd_Comp(inl))<criterium_yPos,1)];
Fin_All_indx_moveEnd_Comp(inl)=[All_indx_moveEnd_Comp(inl)];

elseif ((file_yPos(All_indx_moveEnd_Comp(inl))>criterium_yPos) && (file_yPos(All_indx_moveStart_Comp(inl))<=criterium_yPos))
Fin_All_indx_moveStart_Comp(inl)=[All_indx_moveStart_Comp(inl)];  
Fin_All_indx_moveEnd_Comp(inl)=[All_indx_moveStart_Comp(inl)+find(file_yPos(All_indx_moveStart_Comp(inl):All_indx_moveEnd_Comp(inl))>criterium_yPos,1)-1];    

elseif ((file_yPos(All_indx_moveEnd_Comp(inl))>criterium_yPos) && (file_yPos(All_indx_moveStart_Comp(inl))>criterium_yPos))
try
Fin_All_indx_moveStart_Comp(inl)=[All_indx_moveStart_Comp(inl)+find(file_yPos(All_indx_moveStart_Comp(inl):All_indx_moveEnd_Comp(inl))<criterium_yPos,1)];
Fin_All_indx_moveEnd_Comp(inl)=[All_indx_moveStart_Comp(inl)+find(file_yPos(All_indx_moveStart_Comp(inl):All_indx_moveEnd_Comp(inl))>criterium_yPos,1)-1];    
catch
    try
    Fin_All_indx_moveStart_Comp(inl)=[All_indx_moveStart_Comp(inl)+find(file_yPos(All_indx_moveStart_Comp(inl):All_indx_moveEnd_Comp(inl))<criterium_yPos+6,1)];
    Fin_All_indx_moveEnd_Comp(inl)=[All_indx_moveStart_Comp(inl)+find(file_yPos(All_indx_moveStart_Comp(inl):All_indx_moveEnd_Comp(inl))>criterium_yPos+6,1)-1];    
    catch
        Fin_All_indx_moveStart_Comp(inl)=All_indx_moveStart_Comp(inl);
        Fin_All_indx_moveEnd_Comp(inl)=All_indx_moveEnd_Comp(inl);
    end
    
 end
else Fin_All_indx_moveStart_Comp(inl)=All_indx_moveStart_Comp(inl);%(intersect(find(file_yPos(indx_StimAreaStart)<criterium_yPos), find(file_yPos(indx_StimAreaEnd)<criterium_yPos)))%;%;
Fin_All_indx_moveEnd_Comp(inl)=All_indx_moveEnd_Comp(inl);%(intersect(find(file_yPos(indx_StimAreaStart)<criterium_yPos), find(file_yPos(indx_StimAreaEnd)<criterium_yPos)))%All_indx_moveEnd_Ref;%;

end
end


% try criterium_yPos=mean(file_yPos(intersect(find(file_stimPos~=-1),touch_indx)))+2*std(file_yPos(find(file_stimPos~=-1))); %-4 /-8 ???
% keiner der eckpunkte darf schweben

% Fin_All_indx_moveStart_Ref=All_indx_moveStart_Ref;%(intersect(find(file_yPos(All_indx_moveStart_Ref)<criterium_yPos), find(file_yPos(All_indx_moveEnd_Ref)<criterium_yPos)));
% Fin_All_indx_moveEnd_Ref=All_indx_moveEnd_Ref;%(intersect(find(file_yPos(All_indx_moveStart_Ref)<criterium_yPos), find(file_yPos(All_indx_moveEnd_Ref)<criterium_yPos)));
% Fin_All_indx_moveStart_Comp=All_indx_moveStart_Comp;%(intersect(find(file_yPos(All_indx_moveStart_Comp)<criterium_yPos), find(file_yPos(All_indx_moveEnd_Comp)<criterium_yPos)));
% Fin_All_indx_moveEnd_Comp=All_indx_moveEnd_Comp;%(intersect(find(file_yPos(All_indx_moveStart_Comp)<criterium_yPos), find(file_yPos(All_indx_moveEnd_Comp)<criterium_yPos)));
% 
% fin1=[];
% fin2=[];
% index_fin=[];
%     % die mittleren werte in einem stroke d?fren nicht in der luft liegen
% for f=1:length(Fin_All_indx_moveStart_Comp)
% fin1(f,:)=(quantile([file_yPos(Fin_All_indx_moveStart_Comp(f):Fin_All_indx_moveEnd_Comp(f))],[.90])<criterium_yPos); % 90% der Werte m?ssen unterhalb des Kriterium sein
% fin2(f,:)=(quantile([file_INstimPos(Fin_All_indx_moveStart_Comp(f):Fin_All_indx_moveEnd_Comp(f))],[.90])==0); % 90% der Werte m?ssen 
% index_fin=find((fin1.*fin2)~=0);
% end
% 
% fin_Ref1=[];
% fin_Ref2=[];
% index_fin_Ref=[];
% for f_Ref=1:length(Fin_All_indx_moveStart_Ref)
% fin_Ref1(f_Ref,:)=(quantile([file_yPos(Fin_All_indx_moveStart_Ref(f_Ref):Fin_All_indx_moveEnd_Ref(f_Ref))],[.90])<criterium_yPos);
% fin_Ref2(f_Ref,:)=(quantile([file_INstimPos(Fin_All_indx_moveStart_Ref(f_Ref):Fin_All_indx_moveEnd_Ref(f_Ref))],[.90])==0);
% index_fin_Ref=find((fin_Ref1.*fin_Ref2)~=0);
% end
% 
% Fin_All_indx_moveStart_Comp=Fin_All_indx_moveStart_Comp(index_fin);%All_indx_moveStart_Comp(intersect(find(file_yPos(All_indx_moveStart_Comp)<criterium_yPos), find(file_yPos(All_indx_moveEnd_Comp)<criterium_yPos)));
% Fin_All_indx_moveEnd_Comp=Fin_All_indx_moveEnd_Comp(index_fin);%All_indx_moveEnd_Comp(intersect(find(file_yPos(All_indx_moveStart_Comp)<criterium_yPos), find(file_yPos(All_indx_moveEnd_Comp)<criterium_yPos)));
% Fin_All_indx_moveStart_Ref=Fin_All_indx_moveStart_Ref(index_fin_Ref);%
% Fin_All_indx_moveEnd_Ref=Fin_All_indx_moveEnd_Ref(index_fin_Ref);


%% Bereinigung um zu kurze/ nicht vorhandene Strokes

%1) kein stroke= start und end gl index
Final_All_indx_moveStart_Comp=Fin_All_indx_moveStart_Comp(find(Fin_All_indx_moveStart_Comp~=Fin_All_indx_moveEnd_Comp));
Final_All_indx_moveEnd_Comp=Fin_All_indx_moveEnd_Comp(find(Fin_All_indx_moveStart_Comp~=Fin_All_indx_moveEnd_Comp));

Final_All_indx_moveStart_Ref=Fin_All_indx_moveStart_Ref(find(Fin_All_indx_moveStart_Ref~=Fin_All_indx_moveEnd_Ref));
Final_All_indx_moveEnd_Ref=Fin_All_indx_moveEnd_Ref(find(Fin_All_indx_moveStart_Ref~=Fin_All_indx_moveEnd_Ref));

%2) zu kurzer Stroke= start und ende weniger als 200 msec entfernt

ind_short_Comp=find((file_time(Final_All_indx_moveEnd_Comp)-file_time(Final_All_indx_moveStart_Comp))<=200);
ind_short_Ref=find((file_time(Final_All_indx_moveEnd_Ref)-file_time(Final_All_indx_moveStart_Ref))<=200);

Final_All_indx_moveStart_Comp(ind_short_Comp)=[];
Final_All_indx_moveEnd_Comp(ind_short_Comp)=[];
Final_All_indx_moveStart_Ref(ind_short_Ref)=[];
Final_All_indx_moveEnd_Ref(ind_short_Ref)=[];

% catch
% Final_All_indx_moveStart_Comp=All_indx_moveStart_Comp;
% Final_All_indx_moveEnd_Comp=All_indx_moveEnd_Comp;
% Final_All_indx_moveStart_Ref=All_indx_moveStart_Ref;
% Final_All_indx_moveEnd_Ref=All_indx_moveEnd_Ref;
% 
% ExplError(tCounter)=1;
% end
 %%
%Koordinaten pro Ber?hrungspunkt
kerner=3;%smoothen der Koordinaten
Delta.X.Comp{tCounter}=file_xPos(Final_All_indx_moveEnd_Comp)-file_xPos(Final_All_indx_moveStart_Comp);
Delta.Z.Comp{tCounter}=file_zPos(Final_All_indx_moveEnd_Comp)-file_zPos(Final_All_indx_moveStart_Comp);
Delta.Y.Comp{tCounter}=file_yPos(Final_All_indx_moveEnd_Comp)-file_yPos(Final_All_indx_moveStart_Comp); %hoehe beim phantom
Delta.X.Ref{tCounter}=file_xPos(Final_All_indx_moveEnd_Ref)-file_xPos(Final_All_indx_moveStart_Ref);
Delta.Z.Ref{tCounter}=file_zPos(Final_All_indx_moveEnd_Ref)-file_zPos(Final_All_indx_moveStart_Ref);
Delta.Y.Ref{tCounter}=file_yPos(Final_All_indx_moveEnd_Ref)-file_yPos(Final_All_indx_moveStart_Ref);

M.expl.Comp{tCounter}=-Delta.Z.Comp{tCounter}./Delta.X.Comp{tCounter}; % in phantom_z ist minus vorne und plus hinten, dehalb hier -z/x
M.expl.Ref{tCounter}=-Delta.Z.Ref{tCounter}./Delta.X.Ref{tCounter};

V.expl.Comp{tCounter}=[Delta.X.Comp{tCounter} Delta.Z.Comp{tCounter}];
V.expl.Ref{tCounter}=[Delta.X.Ref{tCounter} Delta.Z.Ref{tCounter}];
% Steigung der Textur
 M.texture(tCounter)=tand(file_TrialOrientation(tCounter)+90);%!!! Ausrichtung ist orthoganal zur Textur
 
% Time for Stroke:
StrokeDuration.Comp{tCounter}=file_time(Final_All_indx_moveEnd_Comp)-file_time(Final_All_indx_moveStart_Comp);
StrokeDuration.Ref{tCounter}=file_time(Final_All_indx_moveEnd_Ref)-file_time(Final_All_indx_moveStart_Ref);

%Explorationswinkel 
 edgesN=[-90:10:90];
Angle.texture{tCounter}=(file_TrialOrientation(tCounter)+90);
Angle.Comp{tCounter}=atand(M.expl.Comp{tCounter});
Angle.hist.Comp(tCounter,:)= histcounts(Angle.Comp{tCounter}(2:end),edgesN);% ab 2tem Stroke
Angle.Ref{tCounter}=atand(M.expl.Ref{tCounter});
Angle.hist.Ref(tCounter,:)= histcounts(Angle.Ref{tCounter}(2:end),edgesN);
 
 % Schnittwinkel
        IntAngle.Comp{tCounter}=atand(abs((M.texture(tCounter)-M.expl.Comp{tCounter})./(1+M.expl.Comp{tCounter}.*M.texture(tCounter)))); %mit abs(?) Vorsicht berechnet immer den kleineren der 2 m?glichen Winkel
        IntAngle.Ref{tCounter}=atand(abs((M.texture(tCounter)-M.expl.Ref{tCounter})./(1+M.expl.Ref{tCounter}.*M.texture(tCounter))));
        IntAngle.Bias.Comp{tCounter}=ones(size(IntAngle.Comp{tCounter}))*90-IntAngle.Comp{tCounter}; % Abwechung von 90?
        IntAngle.Bias.Ref{tCounter}=ones(size(IntAngle.Ref{tCounter}))*90-IntAngle.Ref{tCounter};
        
  edges=[0:10:90];
  IntAngle.Bias.hist.Comp(tCounter,:)= histcounts(IntAngle.Bias.Comp{tCounter}(2:end),edges);
  IntAngle.Bias.hist.Ref(tCounter,:)= histcounts(IntAngle.Bias.Ref{tCounter}(2:end),edges);

  
%% PLOTS        Ist ausgeklammert, da sich zu viele Figuren öffnen
%figure(vp*100+tCounter)
% 
%plot3(file_xPos,file_zPos,file_yPos)
%hold on
% 
% plot3(file_xPos(Final_All_indx_moveStart_Comp),file_zPos(Final_All_indx_moveStart_Comp),file_yPos(Final_All_indx_moveStart_Comp), 'og','MarkerSize',8)
% plot3(file_xPos(Final_All_indx_moveEnd_Comp),file_zPos(Final_All_indx_moveEnd_Comp),file_yPos(Final_All_indx_moveEnd_Comp), 'or')
% plot3(file_xPos(Final_All_indx_moveStart_Ref),file_zPos(Final_All_indx_moveStart_Ref),file_yPos(Final_All_indx_moveStart_Ref), 'pg','MarkerSize',8)
%plot3(file_xPos(Final_All_indx_moveEnd_Ref),file_zPos(Final_All_indx_moveEnd_Ref),file_yPos(Final_All_indx_moveEnd_Ref), 'pr')
% 
% 
% 
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
% title('Positions within 1 Trial')
% 


%FIRST MIDDLE LAST
%first

for i = 1:length(Angle.Comp)
      if isempty(Angle.Comp{1,i})
    Angle.Comp{1,i} = 1;
      end
end
  

for i = 1:length(Angle.Ref)
      if isempty(Angle.Ref{1,i})
    Angle.Ref{1,i} = 1;
      end
end
  

positionindex = find(file_stimPos(:,1)==-1); %alle outerstim postionen rausfiltern. 0 ist links und 1 ist rechts
file_stimPos(positionindex,:) = [];

if file_stimPos(3) == 0
    frequencyfirst = leftStim;
else frequencyfirst = rightStim;
end

if frequencyfirst(1,1) < StimPair(tCounter, 1) | frequencyfirst(1,1) < StimPair(tCounter, 2);
    FirstAngle(tCounter,1) = Angle.Comp{1,tCounter} (1);
else FirstAngle(tCounter,1) = Angle.Ref{1,tCounter} (1);
end


%middle

if file_stimPos(ceil(end/2)) == 0
    frequencymiddle = leftStim;
else frequencymiddle = rightStim;
end


if frequencymiddle(1,1) < StimPair(tCounter, 1) | frequencymiddle(1,1) < StimPair(tCounter, 2)
   MiddleAngle(tCounter,1) = Angle.Comp{1,tCounter} (ceil(end/2));
else MiddleAngle(tCounter,1) = Angle.Ref{1,tCounter} (ceil(end/2));
end

%last

if file_stimPos(end) == 0
    frequencylast = leftStim;
else frequencylast = rightStim;
end

if frequencylast(1,1) < StimPair(tCounter, 1) | frequencylast(1,1) < StimPair(tCounter, 2)
    LastAngle(tCounter,1) = Angle.Comp{1,tCounter} (end);
else LastAngle(tCounter,1) = Angle.Ref{1,tCounter} (end);
end



tCounter=tCounter+1;    %UM DANN BEIM NÄCHSTEN TRIAL WEITERUZMACHEN  ?  
end
end

%% Aaron Results matrix to check

FirstMovement = zeros(360,1);
FlagIndex_Ref = zeros (360,1);
FlagIndex_Comp = zeros (360,1);

for i = 1:360
   if  isempty (IntAngle.Ref{1,i})
       FlagIndex_Ref(i,1) = 1;
   end
end

for i = 1:360
   if  isempty (IntAngle.Comp{1,i})
       FlagIndex_Comp(i,1) = 1;
   end
end



for i = 1:360
    if (StimPair(i,1) == file_refStim(i,1) & FlagIndex_Ref(i,1) == 0)
        FirstMovement (i,1) = IntAngle.Ref{1,i}(1);
    elseif (StimPair(i,1) == file_compStim(i,1) & FlagIndex_Comp(i,1) == 0)
        FirstMovement (i,1) = IntAngle.Comp{1,i}(1);
    end
end


file_Results = zeros(360,5);

file_Results(:,1) = file_trialNr;
file_Results(:,2) = file_TrialQuality;
file_Results(:,3) = FirstMovement;
file_Results(:,4:5) = StimPair;



eval(['Results_VP' num2str(vp) '= file_Results;']);


% DIE DREI von first middle last ZUSAMMENFÜHREN in struct


    move.first = FirstAngle;
    move.middle = MiddleAngle;
    move.last= LastAngle;

%% Save Data
save([savePath, 'S' num2str(vp) '_Stim'], 'StimPair')%,'O_seq','InStimulus')
save([savePath, 'S' num2str(vp) '_timeANDanswer'],'responses','no_cha','time_on_stim1', 'time_on_stim2', 'time_on_Comp', 'time_on_Ref')
save([savePath, 'S' num2str(vp) '_Exploration'],'ExplError', 'M','V', 'Angle','IntAngle', 'Delta', 'StrokeDuration')
save([savePath, 'S' num2str(vp) '_Results'], ['Results_VP' num2str(vp)]);
save([savePath, 'S' num2str(vp) '_Moves'], 'move')
end
