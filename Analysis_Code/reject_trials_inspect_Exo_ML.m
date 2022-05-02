function [EEG,rejtrial] = reject_trials_inspect_Exo_ML(EEG,part_name,exp,rejtrial,i_rej)

i_set = 1; %code changed so now always 1


% /////////////////////////////////////////////////////////////////////////
% `````````````````````````````````````````````````````````````````````````
% ````````````````````````````````````````````````````````````````````````` 
%% Target-aligned epochs
%only when using this specific settings
if strcmpi(exp.setname,'byTargets_ML_v1')  
    
    if strcmpi(part_name,'004') 
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(86)=1;
        EEG.reject.rejthresh(164:165)=1;
        EEG.reject.rejthresh(167)=1;
        EEG.reject.rejthresh(170)=1;
        EEG.reject.rejthresh(179)=1;
        EEG.reject.rejthresh(210)=1;
        EEG.reject.rejthresh(218)=1;
        EEG.reject.rejthresh(243)=1;
        EEG.reject.rejthresh(306)=1;
        EEG.reject.rejthresh(308)=1;
        EEG.reject.rejthresh(310)=1;
        EEG.reject.rejthresh(314:316)=1;
        EEG.reject.rejthresh(318)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
    
        
    elseif strcmpi(part_name,'005') 
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(242)=1;
        EEG.reject.rejthresh(273)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
    
        
    elseif strcmpi(part_name,'006')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(33:34)=1;
        EEG.reject.rejthresh(51)=1;
        EEG.reject.rejthresh(109)=1;
        EEG.reject.rejthresh(149)=1;
        EEG.reject.rejthresh(176)=1;
        EEG.reject.rejthresh(229)=1;
        EEG.reject.rejthresh(245:246)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
    
        
    elseif strcmpi(part_name,'007')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
     
        
    elseif strcmpi(part_name,'008')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(69)=1;
        EEG.reject.rejthresh(241)=1;
        EEG.reject.rejthresh(272)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'009') 
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(311)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
    

    elseif strcmpi(part_name,'010')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(1)=1;
        EEG.reject.rejthresh(53)=1;
        EEG.reject.rejthresh(129)=1;
        EEG.reject.rejthresh(175:176)=1;
        EEG.reject.rejthresh(195)=1;
        EEG.reject.rejthresh(198)=1;
        EEG.reject.rejthresh(297)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'011')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(72)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'012')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(3)=1;
        EEG.reject.rejthresh(64)=1;
        EEG.reject.rejthresh(106)=1;
        EEG.reject.rejthresh(153)=1;
        EEG.reject.rejthresh(181)=1;
        EEG.reject.rejthresh(238)=1;
        EEG.reject.rejthresh(324)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
     
        
    elseif strcmpi(part_name,'013')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(84)=1;
        EEG.reject.rejthresh(267)=1;
        EEG.reject.rejthresh(344)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
     
        
    elseif strcmpi(part_name,'014')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(308)=1;
        EEG.reject.rejthresh(334)=1;
        EEG.reject.rejthresh(355)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
     
        
    elseif strcmpi(part_name,'015')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(156)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
     
        
    elseif strcmpi(part_name,'016')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(169)=1;
        EEG.reject.rejthresh(323)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'017') %need to remove CP1 to keep end of trials
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(154)=1;
        EEG.reject.rejthresh(211)=1;
        EEG.reject.rejthresh(299)=1;
        EEG.reject.rejthresh(301:306)=1;
        EEG.reject.rejthresh(310:338)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
        
        
    elseif strcmpi(part_name,'018')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(147)=1;
        EEG.reject.rejthresh(183)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);   
        
        
    elseif strcmpi(part_name,'019')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(250:251)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);   
        
        
    elseif strcmpi(part_name,'020')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(9)=1;
        EEG.reject.rejthresh(49:50)=1;
        EEG.reject.rejthresh(96)=1;
        EEG.reject.rejthresh(189)=1;
        EEG.reject.rejthresh(261)=1;
        EEG.reject.rejthresh(266)=1;
        EEG.reject.rejthresh(288)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
        
        
    elseif strcmpi(part_name,'021')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(136)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
        
        
    elseif strcmpi(part_name,'022')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(51)=1;
        EEG.reject.rejthresh(114)=1;
        EEG.reject.rejthresh(171)=1;
        EEG.reject.rejthresh(199)=1;
        EEG.reject.rejthresh(203)=1;
        EEG.reject.rejthresh(214)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'023')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);          
     
        
    elseif strcmpi(part_name,'024')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'025')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
     
        
    elseif strcmpi(part_name,'026')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);   
        
        
    elseif strcmpi(part_name,'027')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(95)=1;
        EEG.reject.rejthresh(127:128)=1;
        EEG.reject.rejthresh(242)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);      

        
    elseif strcmpi(part_name,'028')      
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        

    elseif strcmpi(part_name,'029')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(2)=1;
        EEG.reject.rejthresh(123)=1;
        EEG.reject.rejthresh(142)=1;
        EEG.reject.rejthresh(202)=1;
        EEG.reject.rejthresh(228:229)=1;
        EEG.reject.rejthresh(263)=1;
        EEG.reject.rejthresh(266:267)=1;
        EEG.reject.rejthresh(273)=1;
        EEG.reject.rejthresh(277)=1;
        EEG.reject.rejthresh(281)=1;
        EEG.reject.rejthresh(299)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
        

    elseif strcmpi(part_name,'030')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(31)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);         
        

    elseif strcmpi(part_name,'031')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(223)=1;
        EEG.reject.rejthresh(364)=1;
        EEG.reject.rejthresh(368)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);   
        

    elseif strcmpi(part_name,'032')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'033')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        

    elseif strcmpi(part_name,'034')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(155)=1;
        EEG.reject.rejthresh(183)=1;
        EEG.reject.rejthresh(214)=1;
        EEG.reject.rejthresh(281)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'035')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(278)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'036')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(48)=1;
        EEG.reject.rejthresh(138)=1;
        EEG.reject.rejthresh(199)=1;
        EEG.reject.rejthresh(263)=1;
        EEG.reject.rejthresh(292)=1;
        EEG.reject.rejthresh(308)=1;
        EEG.reject.rejthresh(344)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
        
         
    elseif strcmpi(part_name,'037')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(46)=1;
        EEG.reject.rejthresh(82)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
        
        
    elseif strcmpi(part_name,'038')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(132)=1;
        EEG.reject.rejthresh(287)=1;
        EEG.reject.rejthresh(303)=1;
        EEG.reject.rejthresh(321)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
    elseif strcmpi(part_name,'039')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
                  
    elseif strcmpi(part_name,'040')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(13)=1;
        EEG.reject.rejthresh(77)=1;
        EEG.reject.rejthresh(79)=1;
        EEG.reject.rejthresh(93)=1;
        EEG.reject.rejthresh(122)=1;
        EEG.reject.rejthresh(125)=1;
        EEG.reject.rejthresh(131)=1;
        EEG.reject.rejthresh(142)=1;
        EEG.reject.rejthresh(146)=1;
        EEG.reject.rejthresh(151)=1;
        EEG.reject.rejthresh(180)=1;
        EEG.reject.rejthresh(199)=1;
        EEG.reject.rejthresh(205:206)=1;
        EEG.reject.rejthresh(236)=1;
        EEG.reject.rejthresh(249)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
    
     
    end


% /////////////////////////////////////////////////////////////////////////
% `````````````````````````````````````````````````````````````````````````
% ````````````````````````````````````````````````````````````````````````` 
%% Cue-aligned epochs
%only when using this specific settings
elseif strcmpi(exp.setname,'byCues_ML_v1')  
    
     if strcmpi(part_name,'004') 
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(85)=1;
        EEG.reject.rejthresh(126)=1;
        EEG.reject.rejthresh(167)=1;
        EEG.reject.rejthresh(169)=1;
        EEG.reject.rejthresh(172)=1;
        EEG.reject.rejthresh(244)=1;
        EEG.reject.rejthresh(250)=1;
        EEG.reject.rejthresh(297:298)=1;
        EEG.reject.rejthresh(317)=1;
        EEG.reject.rejthresh(321)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
        
    elseif strcmpi(part_name,'005') 
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(243)=1;
        EEG.reject.rejthresh(274)=1;
        EEG.reject.rejthresh(283)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
    
        

     elseif strcmpi(part_name,'006')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(34:35)=1;
        EEG.reject.rejthresh(53)=1;
        EEG.reject.rejthresh(131)=1;
        EEG.reject.rejthresh(139)=1;
        EEG.reject.rejthresh(173)=1;
        EEG.reject.rejthresh(183)=1;
        EEG.reject.rejthresh(186:187)=1;
        EEG.reject.rejthresh(189)=1;
        EEG.reject.rejthresh(198)=1;
        EEG.reject.rejthresh(231:232)=1;
        EEG.reject.rejthresh(246)=1;
        EEG.reject.rejthresh(260)=1;
        EEG.reject.rejthresh(302)=1;
        EEG.reject.rejthresh(317)=1;
        EEG.reject.rejthresh(327:328)=1;
        EEG.reject.rejthresh(339)=1;
        EEG.reject.rejthresh(369)=1;
        EEG.reject.rejthresh(379)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
     
        
     elseif strcmpi(part_name,'007')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
        
        

    elseif strcmpi(part_name,'008')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(69)=1;
        EEG.reject.rejthresh(240)=1;
        EEG.reject.rejthresh(271)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
        
        
        
     elseif strcmpi(part_name,'009') 
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(315)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  

        
        
    elseif strcmpi(part_name,'010')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(1)=1;
        EEG.reject.rejthresh(53)=1;
        EEG.reject.rejthresh(130)=1;
        EEG.reject.rejthresh(176:177)=1;
        EEG.reject.rejthresh(299)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        

     elseif strcmpi(part_name,'011')
        EEG.reject.rejthresh(73)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
       
        
     elseif strcmpi(part_name,'012')
        EEG.reject.rejthresh(65)=1;
        EEG.reject.rejthresh(154)=1;
        EEG.reject.rejthresh(182)=1;
        EEG.reject.rejthresh(239)=1;
        EEG.reject.rejthresh(326)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
     elseif strcmpi(part_name,'013')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(170)=1;
        EEG.reject.rejthresh(349)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
     
        
     elseif strcmpi(part_name,'014')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(334)=1;
        EEG.reject.rejthresh(355)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
     elseif strcmpi(part_name,'015')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);   
       
     
     elseif strcmpi(part_name,'016')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(165)=1;
        EEG.reject.rejthresh(319)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
           
     
     elseif strcmpi(part_name,'017') %need to remove CP1 to keep end of trials
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(212)=1;
        EEG.reject.rejthresh(298)=1;
        EEG.reject.rejthresh(301)=1;
        EEG.reject.rejthresh(304:309)=1;
        EEG.reject.rejthresh(313:341)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
  
        
     elseif strcmpi(part_name,'018')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(184)=1;
        EEG.reject.rejthresh(231)=1;
        EEG.reject.rejthresh(233)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
        
        
     elseif strcmpi(part_name,'019')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(253:254)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
    
     elseif strcmpi(part_name,'021')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(34)=1;
        EEG.reject.rejthresh(136)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
     
        
     elseif strcmpi(part_name,'022')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(52)=1;
        EEG.reject.rejthresh(115)=1;
        EEG.reject.rejthresh(172)=1;
        EEG.reject.rejthresh(202)=1;
        EEG.reject.rejthresh(247)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
     elseif strcmpi(part_name,'023')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);         
     
        
     elseif strcmpi(part_name,'024')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(318)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
        
        
     elseif strcmpi(part_name,'025')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
     
        
     elseif strcmpi(part_name,'026')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
        
        
      elseif strcmpi(part_name,'027')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(95)=1;
        EEG.reject.rejthresh(127)=1;
        EEG.reject.rejthresh(317)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);      
        
        
     elseif strcmpi(part_name,'028')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(69)=1;
        EEG.reject.rejthresh(79)=1;
        EEG.reject.rejthresh(144)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);         
        
        
     elseif strcmpi(part_name,'029')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(2)=1;
        EEG.reject.rejthresh(57:58)=1;
        EEG.reject.rejthresh(64)=1;
        EEG.reject.rejthresh(77)=1;
        EEG.reject.rejthresh(88)=1;
        EEG.reject.rejthresh(92)=1;
        EEG.reject.rejthresh(109)=1;
        EEG.reject.rejthresh(123)=1;
        EEG.reject.rejthresh(142)=1;
        EEG.reject.rejthresh(193)=1;
        EEG.reject.rejthresh(203)=1;
        EEG.reject.rejthresh(222)=1;
        EEG.reject.rejthresh(230)=1;
        EEG.reject.rejthresh(231)=1;
        EEG.reject.rejthresh(261)=1;
        EEG.reject.rejthresh(265:266)=1;
        EEG.reject.rejthresh(269:270)=1;
        EEG.reject.rejthresh(277)=1;
        EEG.reject.rejthresh(281)=1;
        EEG.reject.rejthresh(285)=1;
        EEG.reject.rejthresh(303)=1;
        EEG.reject.rejthresh(312)=1;
        EEG.reject.rejthresh(314)=1;
        EEG.reject.rejthresh(322)=1;
        EEG.reject.rejthresh(324)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
        
      
     elseif strcmpi(part_name,'030')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(15)=1;
        EEG.reject.rejthresh(49:50)=1;
        EEG.reject.rejthresh(56)=1;
        EEG.reject.rejthresh(58:59)=1;
        EEG.reject.rejthresh(65)=1;
        EEG.reject.rejthresh(74)=1;
        EEG.reject.rejthresh(101)=1;
        EEG.reject.rejthresh(110)=1;
        EEG.reject.rejthresh(137)=1;
        EEG.reject.rejthresh(140)=1;
        EEG.reject.rejthresh(150)=1;
        EEG.reject.rejthresh(158)=1;
        EEG.reject.rejthresh(191)=1;
        EEG.reject.rejthresh(221)=1;
        EEG.reject.rejthresh(225)=1;
        EEG.reject.rejthresh(241)=1;
        EEG.reject.rejthresh(252)=1;
        EEG.reject.rejthresh(266)=1;
        EEG.reject.rejthresh(286:287)=1;
        EEG.reject.rejthresh(316)=1;
        EEG.reject.rejthresh(327)=1;
        EEG.reject.rejthresh(332)=1;
        EEG.reject.rejthresh(328)=1;
        EEG.reject.rejthresh(337)=1;
        EEG.reject.rejthresh(340)=1;
        EEG.reject.rejthresh(355)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);         
        
        
     elseif strcmpi(part_name,'031')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(223)=1;
        EEG.reject.rejthresh(341)=1;
        EEG.reject.rejthresh(364)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);   

        
      elseif strcmpi(part_name,'032')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
        
        
     elseif strcmpi(part_name,'033')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(25)=1;
        EEG.reject.rejthresh(145)=1;
        EEG.reject.rejthresh(207)=1;
        EEG.reject.rejthresh(317)=1;
        EEG.reject.rejthresh(325)=1;
        EEG.reject.rejthresh(353)=1;
        EEG.reject.rejthresh(363)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 

        
     elseif strcmpi(part_name,'034')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(159)=1;
        EEG.reject.rejthresh(188)=1;
        EEG.reject.rejthresh(219)=1;
        EEG.reject.rejthresh(286)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
     elseif strcmpi(part_name,'035')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(288)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
      
     elseif strcmpi(part_name,'036')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(48)=1;
        EEG.reject.rejthresh(138)=1;
        EEG.reject.rejthresh(190)=1;
        EEG.reject.rejthresh(200)=1;
        EEG.reject.rejthresh(212)=1;
        EEG.reject.rejthresh(264)=1;
        EEG.reject.rejthresh(293)=1;
        EEG.reject.rejthresh(309)=1;
        EEG.reject.rejthresh(345)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
        
         
     elseif strcmpi(part_name,'037')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(46)=1;
        EEG.reject.rejthresh(82)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
             
        
     elseif strcmpi(part_name,'038')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(178)=1;
        EEG.reject.rejthresh(292)=1;
        EEG.reject.rejthresh(325)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
        
        
     elseif strcmpi(part_name,'039')
%         EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        rejtrial(i_set,i_rej).ids = [];
%         EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);   
                

     elseif strcmpi(part_name,'040')
        EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
        EEG.reject.rejthresh(1)=1;
        EEG.reject.rejthresh(13)=1;
        EEG.reject.rejthresh(64)=1;
        EEG.reject.rejthresh(73)=1;
        EEG.reject.rejthresh(75)=1;
        EEG.reject.rejthresh(77)=1;
        EEG.reject.rejthresh(91)=1;
        EEG.reject.rejthresh(119)=1;
        EEG.reject.rejthresh(122)=1;
        EEG.reject.rejthresh(128)=1;
        EEG.reject.rejthresh(139)=1;
        EEG.reject.rejthresh(143)=1;
        EEG.reject.rejthresh(148)=1;
        EEG.reject.rejthresh(177)=1;
        EEG.reject.rejthresh(196)=1;
        EEG.reject.rejthresh(202:203)=1;
        EEG.reject.rejthresh(233)=1;
        EEG.reject.rejthresh(246)=1;
        rejtrial(i_set,i_rej).ids = find(EEG.reject.rejthresh==1);
        EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
        
        
     end
     

% ````````````````````````````````````````````````````````````````````````````````````````````
% ````````````````````````````````````````````````````````````````````````````````````````````
end

% ````````````````````````````````````````````````````````````````````````````````````````````
% ````````````````````````````````````````````````````````````````````````````````````````````


% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% save rejected trials
EEG.rejtrial = rejtrial;
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::














