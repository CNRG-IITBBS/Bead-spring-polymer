%%%%%%%%%%%%%%%%%%% Input Parameters wt%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%
clear all; clc; close all;

r=0.02; % addition of t.f
g=0.25; % remove t.f 
e_no=44;            % number of enhancers
%% Import data from text file

opts = delimitedTextImportOptions("NumVariables", 16);
opts.DataLines = [2, Inf];
opts.Delimiter = " ";

opts.VariableNames = ["cl", "r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "r10", "r11", "r12", "r13"];
opts.VariableTypes = ["double","double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
ratesbindingfinale = readtable("rates_binding_finale.csv", opts);
ratesunbindingfinale = readtable("rates_unbinding_finale.csv", opts);
clear opts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for sox9%%%%%%%%%%%%%%%
ep21=table2array(ratesbindingfinale(:,2));
ep22=table2array(ratesbindingfinale(:,3));
ep23=table2array(ratesbindingfinale(:,4));
ep24=table2array(ratesbindingfinale(:,5));
ep25=table2array(ratesbindingfinale(:,6));
ep26=table2array(ratesbindingfinale(:,7));
ep27=table2array(ratesbindingfinale(:,8));
ep28=table2array(ratesbindingfinale(:,9));
ep29=table2array(ratesbindingfinale(:,10));
ep210=table2array(ratesbindingfinale(:,11));
ep211=table2array(ratesbindingfinale(:,12));
ep212=table2array(ratesbindingfinale(:,13));
ep213=table2array(ratesbindingfinale(:,14));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epub21=table2array(ratesunbindingfinale(:,2));
epub22=table2array(ratesunbindingfinale(:,3));
epub23=table2array(ratesunbindingfinale(:,4));
epub24=table2array(ratesunbindingfinale(:,5));
epub25=table2array(ratesunbindingfinale(:,6));
epub26=table2array(ratesunbindingfinale(:,7));
epub27=table2array(ratesunbindingfinale(:,8));
epub28=table2array(ratesunbindingfinale(:,9));
epub29=table2array(ratesunbindingfinale(:,10));
epub210=table2array(ratesunbindingfinale(:,11));
epub211=table2array(ratesunbindingfinale(:,12));
epub212=table2array(ratesunbindingfinale(:,13));
epub213=table2array(ratesunbindingfinale(:,14));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%for kcnj2%%%%%%%%%%%%%%%%

ep1=zeros(1,e_no);
epub1=zeros(1,e_no);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epp1={ep1,ep1,ep1,ep1,ep1,ep1,ep1,ep1,ep1,ep1,ep1,ep1,ep1}; 
eppub1={epub1,epub1,epub1,epub1,epub1,epub1,epub1,epub1,epub1,epub1,epub1,epub1,epub1};
epp2={ep21,ep22,ep23,ep24,ep25,ep26,ep27,ep28,ep29,ep210,ep211,ep212,ep213};
eppub2={epub21,epub22,epub23,epub24,epub25,epub26,epub27,epub28,epub29,epub210,epub211,epub212,epub213};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sites=1;

runs=10000; % total simulation runs
fileID1=fopen('trans_km_scaling_case_wt_1.5sigma_rgratio_5_together_off.csv','w');
fprintf(fileID1,'%s %s %s %s %s\n','netpolp1','netpolp2','netenp1','netenp2');
start=1000;
snap=linspace(start,runs,1000);

%%%%%%%%%%%%%%%%%%%start the simulation %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%start the simulation %%%%%%%%%%%
for nj=1:1
    cmax=44;
    cmax2=9;
    nj;
    kout1=0;
    kout2=0;
    P1=zeros(1,sites);
    P2=zeros(1,sites);
    p1=zeros(1,e_no);
    p2=zeros(1,e_no);
    tote1=0;
    tote2=0;
    adpol1=0;adpol2=0; rempol1=0; rempol2=0; adpolep1=0;adpolep2=0;time=0;
    pp=0;
    vv='reading_wt'+string(cmax)+string(nj)+'.csv'
    fileID4=fopen(vv,'w');
    for i=1:runs
        %%%%%%%%%%%%%%%%
        x=kout1+1;
        y=kout2+1;
        x=sum(p1)+1;
        y=sum(p2)+1;
        ep1=cell2mat(epp1(x));
        epub1=cell2mat(eppub1(x));
        ep2=cell2mat(epp2(y));
        epub2=cell2mat(eppub2(y));
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          
        


        i;
        m=0;
        n=0;
        rate=0;
        mm=0;
        st_rt={};
        

       some1=sum(p1);
       some2=sum(p2);
       some3=P1(1);
       some4=P2(1);
        %%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        if some3<cmax2
            if P1(1)==0 
                st_rt{1}=r;
                st_rt{2}=0;
           
        
            else
                st_rt{1}=r;
                st_rt{2}=g;
           
            end
        elseif some3>=cmax2
            st_rt{1}=0;
            st_rt{2}=g;
        end 

        if some4<cmax2

            if P2(1)==0
                st_rt{3}=r;
                st_rt{4}=0;
           
            else
                st_rt{3}=r;
                st_rt{4}=g;
            end    
        elseif some4>=cmax2
            st_rt{3}=0;
            st_rt{4}=g;
        end
      
        fact1=1.0;
        fact2=1.0;
        fact3=1.0;
        fact4=1.0;

        for j=1:e_no % number of enhancers
            m=m+1;
            if p1(j)==1 && p2(j)==0
                
                rate=rate+epub1(j)*fact3; % where enhancers are present sum them

                st_rt{5}(m)=epub1(j)*fact3; % unbinding rates of thes enhancers
            else

                st_rt{5}(m)=0;
            end
            n=n+1;
            if p2(j)==1 && p1(j)==0
                fact4=1;%(1.037*sum(p2)-0.336);
                rate=rate+epub2(j)*fact4; % where enhancers are present
                st_rt{6}(n)=epub2(j)*fact4; % unbinding rates of the enhancers present at p2
            else
                st_rt{6}(n)=0;

            end
            mm=mm+1;
            if p1(j)==0 && p2(j)==0
                if some1<cmax
                   
                    rate=rate+ep1(j)*fact1; % add to the total rate
                    st_rt{7}(mm)=ep1(j)*fact1; %here i am binding rates at p1
                elseif some1>=cmax
                    st_rt{7}(mm)=0;

                end
                if some2<cmax
                    %fact2=1/sum(ep2);
                    rate=rate+ep2(j)*fact2;
                    st_rt{8}(mm)=ep2(j)*fact2; % binding rates at p2
                elseif some2>=cmax
                    st_rt{8}(mm)=0;
                end
            else
                st_rt{7}(mm)=0;
                st_rt{8}(mm)=0;

            end
        end

        R=st_rt{1}+st_rt{2}+st_rt{3}+st_rt{4}+rate; %  sum rates of all the reactions
        r2=R*rand(); % generate random no b/w 0 and R
        r1=rand(); % this one for time
        t=(1/R)*log(1/r1); % time form exponential distributio
        time=time+t; % time at which next event happen



        % creating an cell array
        loc=randi(sites); % random site location for adding pol
        rate_total=st_rt{1}+st_rt{2}+st_rt{3}+st_rt{4};

        % snapshots%%%%%%%%%%%%%%%%%%%%%

        if i==10000+1*pp
            pp=pp+1;

            fprintf(fileID4,'%5.3f %5.3f %5.3f %5.3f\n',P1(1),P2(1),sum(p1),sum(p2));
        end




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if  (0<r2) && r2<=st_rt{1} && P1(loc)<cmax2
            "step1";
            P1(loc)=P1(loc)+1;
            adpol1=adpol1+1;

        elseif (st_rt{1}<r2) && (r2<=(st_rt{1}+st_rt{2})) && P1(loc)>0
            P1(loc)=P1(loc)-1; % romove pol at the site loc
            rempol1=rempol1+1;
            "step2";
        elseif ((st_rt{1}+st_rt{2})<r2) && (r2<=(st_rt{1}+st_rt{2}+st_rt{3})) &&  P2(loc)<cmax2
            P2(loc)=P2(loc)+1;
            adpol2=adpol2+1;
            "step3";
        elseif ((st_rt{1}+st_rt{2}+st_rt{3})<r2) &&  (r2<=(rate_total)) && P2(loc)>0
            P2(loc)=P2(loc)-1;
            rempol2=rempol2+1;
            "step4";
        else  % you are here
            "step5";
            % t3
            pos={};
            for im=1:4
                if length(st_rt{im+4})>0
                    for j=1:length(st_rt{im+4})

                        if (rate_total<r2)&& (r2<=(rate_total+st_rt{im+4}(j)))
                            pos{im}(j)=1;

                        else
                            pos{im}(j)=0;
                        end
                        rate_total=rate_total+st_rt{im+4}(j);


                    end

                elseif length(st_rt{im+4})==0
                    for j=1:e_no
                        pos{im}(j)=0;
                    end


                end
            end

            for ii= 1:length(pos)
                if (rem(ii,2)==1) && (sum(pos{ii})>0)
                    for jj=1:length(pos{ii})
                        if ii==1 && pos{ii}(jj)==1  % ubinding of enhancer from p1
                            p1(jj)=p1(jj)-1; % remove enhancer at jj position
                            kout1=kout1-1;

                        elseif ii==3 && pos{ii}(jj)==1 && P1(1)<cmax2
                            p1(jj)=p1(jj)+1; % add enhancer at jj position
                            %  P1(loc)=P1(loc)+1; % add pol at loc
                            kout1=kout1+1;
                            % adpolep1=adpolep1+1;

                        end
                    end

                elseif (rem(ii,2)==0) && (sum(pos{ii})>0)
                    for jj=1:length(pos{ii})
                        if ii==2 && pos{ii}(jj)==1 % ubinding of enhancer from p1

                            p2(jj)=p2(jj)-1; % remove enhancer at jj position
                            kout2=kout2-1;

                        elseif ii==4 && pos{ii}(jj)==1 && P2(1)<cmax2
                            p2(jj)=p2(jj)+1; % add enhancer at jj position
                            % P2(loc)=P2(loc)+1; % add pol at loc
                            kout2=kout2+1;
                            %  adpolep2=adpolep2+1;
                        end
                    end
                end
            end


        end


        tote1=tote1+sum(p1);
        tote2=tote2+sum(p2);



    end

   netpolp1=adpol1-rempol1;
   netpolp2=adpol2-rempol2;
    fprintf(fileID1,'%5.3f %5.3f %5.3f %5.3f \n',netpolp1,netpolp2,kout1,kout2);
     
end
fclose(fileID1);


                                                                                                                                                                                                                                                                                                                     
