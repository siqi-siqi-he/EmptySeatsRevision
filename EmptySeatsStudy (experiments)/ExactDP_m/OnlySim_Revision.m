for noExtra = [0,1] %=0: selling product 3, =1: not selling product 3 = the price of product 3 = inf
    for choice = [0]
        if noExtra == 0
            i = 2;
        else
            i = 30;
        end
        if choice == 0
            sheet = 16;
            folder = 'SensitivityAnalysisResults/HOMOG_Revision/';
        else
            error('choice parameter not feasible');
        end
        for c = [8] % Never chose an odd number!
            T = 2*c;
            for a0 = [-1,-0.5,0,0.5,1]
                nsim = 100;
                file_name = sprintf('%sDP_par_simulation%d_choice%d_capacity%d_nopurch%d_noExtra%d_a3%d_Revision_median.mat',folder,nsim,choice,c,a0,noExtra,a3(1));
                load(file_name)

                % Get random decision variables
                rng(987654);
                nsim = 100;
                dec = sample_dec(T+1,nsim);
                %use for loop for loading CRNs
                for n =1:nsim
                    path = sprintf('C:/Users/davinah/Dropbox/Empty_Seats/handover_code/revision/Fig8_values/RandomNumbers/randomnumbers_sim%d.txt',n-1);
                    col = load(path);
                    dec(2:end,n) = col; 
                end
                tstart = tic;
                [rev_nt,policy,state,tickets_n,seats2_n,seats3_n] = simulation(nsim,dec,T,c,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,V,noExtra);
                t_sim = toc(tstart);
                stderr = std(sum(rev_nt,2));
                rev = mean(sum(rev_nt,2))

                %save data
                file_name = sprintf('%sDP_par_simulation%d_choice%d_capacity%d_nopurch%d_noExtra%d_a3%d_Revision_median.mat',folder,nsim,choice,c,a0,noExtra,a3(1));
                save(file_name,'V','Vroot','t_training','a0','a1','a2','a3','b1','b2','b3','tau1','tau2','tau3','c','T','rev_nt','rev','policy','tickets_n','seats2_n','seats3_n');
                filename = 'Results_DP_Revision.xlsx';
%                 writematrix('',filename,'Sheet',sheet,'Range','A'+string(i));
%                 writematrix(c,filename,'Sheet',sheet,'Range','B'+string(i));
%                 writematrix(Vroot,filename,'Sheet',sheet,'Range','C'+string(i));
%                 writematrix(t_training,filename,'Sheet',sheet,'Range','D'+string(i));
                writematrix(rev,filename,'Sheet',sheet,'Range','E'+string(i));
                writematrix(t_sim,filename,'Sheet',sheet,'Range','F'+string(i));
                writematrix(nsim,filename,'Sheet',sheet,'Range','G'+string(i));
%                 writematrix(a0,filename,'Sheet',sheet,'Range','H'+string(i));
%                 writematrix(0,filename,'Sheet',sheet,'Range','I'+string(i));
                writematrix(stderr,filename,'Sheet',sheet,'Range','J'+string(i));
                writematrix(mean(tickets_n),filename,'Sheet',sheet,'Range','K'+string(i)+':N'+string(i));
                writematrix(mean(seats2_n),filename,'Sheet',sheet+1,'Range',string(i)+':'+string(i));
                writematrix(mean(seats3_n),filename,'Sheet',sheet+2,'Range',string(i)+':'+string(i));
                i = i+1;
            end
        end
    end
end
