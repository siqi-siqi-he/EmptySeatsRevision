%% DemandCapAnalysis %%
for noExtra = [0,1] %=0: selling product 3, =1: not selling product 3 = the price of product 3 = inf
    for choice = [0]
        if noExtra == 0
            i = 2;
        else
            i = 7;
        end
        if choice == 0
            sheet = 1;
        elseif choice == 1
            sheet = 4;
        elseif choice == 2
            sheet = 7;
        end
        for c = [8] % Never chose an odd number!
            T = 2*c;
            % Generate values in MATLAB
            if choice == 0
                [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_homo_seats(c);
                folder = 'SensitivityAnalysisResults/DemandCap/HOMOG/';
            elseif choice == 1
                [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_aisle_window_seats(c);
                folder = 'SensitivityAnalysisResults/DemandCap/AW/';
            elseif choice == 2
                [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_increasing_seats(c);
                folder = 'SensitivityAnalysisResults/DemandCap/HET/';
            end
            disp([a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3]);

            for a0 = [-1,-0.5,0,0.5,1]
                % Initialize V table with zeros: V[seat arrangement,reservation,time]
                V = zeros(2^c + 1,c + 2,T + 1);
                SLF = zeros(2^c + 1,c + 2,T + 1);
                SLF(:,1,:) = 0;
                V(:,1,:) = -9999;
                V(1,:,:) = -9999;
                SLFhelp = zeros(2^c + 1,c + 2);
                Vhelp = zeros(2^c + 1,c + 2);
                Vhelp(:,1) = -9999;
                Vhelp(1,:) = -9999;
                policy = zeros(2^c + 1,c + 2,1+2*c,T+1);
                policyhelp = zeros(2^c + 1,c+2,1+2*c);

                % Learning process
                Y = sum(c)+2;
                tstart = tic;
                for t = 2:T+1
                    disp(t-1);
                    parfor D = 3:(2^c+1)
                        x = DtoB0(D,c);
                        for y = 3:Y %Y(D) does not work due to parfor, thus another if clause
                            if y>sum(x)+2
                                Vhelp(D,y) = -9999;
                            else
                                % Optimize prices on expectation value
                                [ofv,r1,r2j,r3j,SLFx] = opt_DP(V,c,x,y,t,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,noExtra,SLF);
                                Vhelp(D,y) = ofv;
                                SLFhelp(D,y) = SLFx;
                                policyhelp(D,y,:) = [r1;r2j;r3j];
                            end
                        end
                    end
                    policy(:,:,:,t) = policyhelp;
                    V(:,:,t) = Vhelp;
                    SLF(:,:,t) = SLFhelp;
                end
                t_training = toc(tstart);

                x = ones(c,1);
                y = sum(x) + 2;
                t = T+1;
                disp(V(BtoD0(x),y,t));
                Vroot = V(BtoD0(x),y,t);
                SLFoverall = SLF(BtoD0(x),y,t);
                disp(SLFoverall/c);
                policy2 = policy;

                % Get random decision variables
                rng(987654);
                nsim = 100;
                dec = sample_dec(T+1,nsim);
                tstart = tic;
                [rev_nt,policy,state,tickets_n,seats2_n,seats3_n] = simulation(nsim,dec,T,c,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,V,noExtra);
                t_sim = toc(tstart);
                stderr = std(sum(rev_nt,2));
                rev = mean(sum(rev_nt,2))

                %save data
                file_name = sprintf('%sDP_par_simulation%d_choice%d_capacity%d_nopurch%d_noExtra%d.mat',folder,nsim,choice,c,a0,noExtra);
                save(file_name,'V','Vroot','t_training','SLF','policy2','a0','a1','a2','a3','b1','b2','b3','tau1','tau2','tau3','c','T','rev_nt','rev','policy','tickets_n','seats2_n','seats3_n');
                filename = 'Results_DP.xlsx';
                writematrix(noExtra,filename,'Sheet',sheet,'Range','A'+string(i));
                writematrix(c,filename,'Sheet',sheet,'Range','B'+string(i));
                writematrix(Vroot,filename,'Sheet',sheet,'Range','C'+string(i));
                writematrix(t_training,filename,'Sheet',sheet,'Range','D'+string(i));
                writematrix(rev,filename,'Sheet',sheet,'Range','E'+string(i));
                writematrix(t_sim,filename,'Sheet',sheet,'Range','F'+string(i));
                writematrix(nsim,filename,'Sheet',sheet,'Range','G'+string(i));
                writematrix(a0,filename,'Sheet',sheet,'Range','H'+string(i));
                writematrix(SLFoverall/c,filename,'Sheet',sheet,'Range','I'+string(i));
                writematrix(stderr,filename,'Sheet',sheet,'Range','J'+string(i));
                writematrix(mean(tickets_n),filename,'Sheet',sheet,'Range','K'+string(i)+':N'+string(i));
                writematrix(mean(seats2_n),filename,'Sheet',sheet+1,'Range',string(i)+':'+string(i));
                writematrix(mean(seats3_n),filename,'Sheet',sheet+2,'Range',string(i)+':'+string(i));
                i = i+1;
            end
        end
    end
end

% %% NestDiss %%
% for contrary = [0,1]
% for noExtra = [1] %=0: selling product 3, =1: not selling product 3 = the price of product 3 = inf
%     for choice = [1,2]
%         if noExtra == 0
%             if contrary == 0
%                 i = 76;
%             else
%                 i = 38;
%             end
%         else
%             if contrary == 0
%                 i = 82;
%             else
%                 i = 47;
%             end
%         end
%         if choice == 0
%             sheet = 1;
%         elseif choice == 1
%             sheet = 4;
%         elseif choice == 2
%             sheet = 7;
%         end
%         for c = [8] % Never chose an odd number!
%             T = 2*c;
% 
%             % Generate values in MATLAB
%             if choice == 0
%                 [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_homo_seats(c);
%                 folder = 'SensitivityAnalysisResults/LevelDissim23/HOMOG/';
%             elseif choice == 1
%                 [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_aisle_window_seats(c);
%                 folder = 'SensitivityAnalysisResults/LevelDissim23/AW/';
%             elseif choice == 2
%                 [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_increasing_seats(c);
%                 folder = 'SensitivityAnalysisResults/LevelDissim23/HET/';
%             end
%             disp([a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3]);
%             if contrary == 1
%                 taus = [[0.1;0.9],[0.2;0.8],[0.3;0.7],[0.4;0.6],[0.5;0.5],[0.6;0.4],[0.7;0.3],[0.8;0.2],[0.9;0.1]];
%             else 
%                 taus = [0.25*[tau2;tau3],0.5.*[tau2;tau3],0.75.*[tau2;tau3],1.*[tau2;tau3],1.25.*[tau2;tau3],1.5.*[tau2;tau3]];
%             end
% 
%             for tau = taus
%                 tau2 = tau(1);
%                 tau3 = tau(2);
% 
%                 % Initialize V table with zeros: V[seat arrangement,reservation,time]
%                 V = zeros(2^c + 1,c + 2,T + 1);
%                 V(:,1,:) = -9999;
%                 V(1,:,:) = -9999;
%                 Vhelp = zeros(2^c + 1,c + 2);
%                 Vhelp(:,1) = -9999;
%                 Vhelp(1,:) = -9999;
%                 policy = zeros(2^c + 1,c + 2,1+2*c,T+1);
%                 policyhelp = zeros(2^c + 1,c+2,1+2*c);
% 
%                 % Learning process
%                 Y = sum(c)+2;
%                 tstart = tic;
%                 for t = 2:T+1
%                     disp(t-1);
%                     parfor D = 3:(2^c+1)
%                         x = DtoB0(D,c);
%                         for y = 3:Y 
%                             if y>sum(x)+2
%                                 Vhelp(D,y) = -9999;
%                             else
%                                 [ofv,r1,r2j,r3j] = opt_DP(V,c,x,y,t,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,noExtra);
%                                 Vhelp(D,y) = ofv;
%                                 policyhelp(D,y,:) = [r1;r2j;r3j];
%                             end
%                         end
%                     end
%                     policy(:,:,:,t) = policyhelp;
%                     V(:,:,t) = Vhelp;
%                 end
%                 t_training = toc(tstart);
% 
%                 x = ones(c,1);
%                 y = sum(x) + 2;
%                 t = T+1;
%                 disp(V(BtoD0(x),y,t));
%                 Vroot = V(BtoD0(x),y,t);
%                 policy2 = policy;
% 
%                 % Get random decision variables
%                 rng(987654);
%                 nsim = 100;
%                 dec = sample_dec(T+1,nsim);
%                 tstart = tic;
%                 [rev_nt,policy,state,tickets_n,seats2_n,seats3_n] = simulation(nsim,dec,T,c,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,V,noExtra);
%                 t_sim = toc(tstart);
%                 stderr = std(sum(rev_nt,2));
%                 rev = mean(sum(rev_nt,2))
% 
%                 %save data
%                 file_name = sprintf('%sDP_par_simulation%d_choice%d_capacity%d_nopurch%d_nestdiss%d_%d_noExtra%d.mat',folder,nsim,choice,c,a0,tau2,tau3,noExtra);
%                 save(file_name,'V','Vroot','t_training','policy2','a0','a1','a2','a3','b1','b2','b3','tau1','tau2','tau3','c','T','rev_nt','rev','policy','tickets_n','seats2_n','seats3_n');
%                 filename = 'Results_DP.xlsx';
%                 writematrix('',filename,'Sheet',sheet,'Range','A'+string(i));
%                 writematrix(c,filename,'Sheet',sheet,'Range','B'+string(i));
%                 writematrix(Vroot,filename,'Sheet',sheet,'Range','C'+string(i));
%                 writematrix(t_training,filename,'Sheet',sheet,'Range','D'+string(i));
%                 writematrix(rev,filename,'Sheet',sheet,'Range','E'+string(i));
%                 writematrix(t_sim,filename,'Sheet',sheet,'Range','F'+string(i));
%                 writematrix(nsim,filename,'Sheet',sheet,'Range','G'+string(i));
%                 writematrix(a0,filename,'Sheet',sheet,'Range','H'+string(i));
% %                 writematrix(tau1,filename,'Sheet',sheet,'Range','I'+string(i));
%                 writematrix(stderr,filename,'Sheet',sheet,'Range','J'+string(i));
%                 writematrix(mean(tickets_n),filename,'Sheet',sheet,'Range','K'+string(i)+':N'+string(i));
%                 writematrix(mean(seats2_n),filename,'Sheet',sheet+1,'Range',string(i)+':'+string(i));
%                 writematrix(mean(seats3_n),filename,'Sheet',sheet+2,'Range',string(i)+':'+string(i));
%                 i = i+1;
%             end
%         end
%     end
% end
% end
% 
% % %% Differences of a1 %%
% % for noExtra = [1] %=0: selling product 3, =1: not selling product 3 = the price of product 3 = inf
% %     for choice = [0] %,1,2
% %         if noExtra == 0
% %             i = 58;
% %         else
% %             i = 66;
% %         end
% %         if choice == 0
% %             sheet = 1;
% %         elseif choice == 1
% %             sheet = 4;
% %         elseif choice == 2
% %             sheet = 7;
% %         end
% %         for c = [8] % Never chose an odd number!
% %             T = 2*c;
% % 
% %             % Generate values in MATLAB
% %             if choice == 0
% %                 [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_homo_seats(c);
% %                 folder = 'SensitivityAnalysisResults/LevelDissim23/HOMOG/';
% %             elseif choice == 1
% %                 [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_aisle_window_seats(c);
% %                 folder = 'SensitivityAnalysisResults/LevelDissim23/AW/';
% %             elseif choice == 2
% %                 [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_increasing_seats(c);
% %                 folder = 'SensitivityAnalysisResults/LevelDissim23/HET/';
% %             end
% %             disp([a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3]);
% %             a1s = [0.05:0.05:0.4];
% % 
% %             for a1 = a1s
% % 
% %                 % Initialize V table with zeros: V[seat arrangement,reservation,time]
% %                 V = zeros(2^c + 1,c + 2,T + 1);
% %                 V(:,1,:) = -9999;
% %                 V(1,:,:) = -9999;
% %                 Vhelp = zeros(2^c + 1,c + 2);
% %                 Vhelp(:,1) = -9999;
% %                 Vhelp(1,:) = -9999;
% %                 policy = zeros(2^c + 1,c + 2,1+2*c,T+1);
% %                 policyhelp = zeros(2^c + 1,c+2,1+2*c);
% % 
% %                 % Learning process
% %                 Y = sum(c)+2;
% %                 tstart = tic;
% %                 for t = 2:T+1
% %                     disp(t-1);
% %                     parfor D = 3:(2^c+1)
% %                         x = DtoB0(D,c);
% %                         for y = 3:Y 
% %                             if y>sum(x)+2
% %                                 Vhelp(D,y) = -9999;
% %                             else
% %                                 [ofv,r1,r2j,r3j] = opt_DP(V,c,x,y,t,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,noExtra);
% %                                 Vhelp(D,y) = ofv;
% %                                 policyhelp(D,y,:) = [r1;r2j;r3j];
% %                             end
% %                         end
% %                     end
% %                     policy(:,:,:,t) = policyhelp;
% %                     V(:,:,t) = Vhelp;
% %                 end
% %                 t_training = toc(tstart);
% % 
% %                 x = ones(c,1);
% %                 y = sum(x) + 2;
% %                 t = T+1;
% %                 disp(V(BtoD0(x),y,t));
% %                 Vroot = V(BtoD0(x),y,t);
% %                 policy2 = policy;
% % 
% %                 % Get random decision variables
% %                 rng(987654);
% %                 nsim = 100;
% %                 dec = sample_dec(T+1,nsim);
% %                 tstart = tic;
% %                 [rev_nt,policy,state,tickets_n,seats2_n,seats3_n] = simulation(nsim,dec,T,c,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,V,noExtra);
% %                 t_sim = toc(tstart);
% %                 stderr = std(sum(rev_nt,2));
% %                 rev = mean(sum(rev_nt,2))
% % 
% %                 %save data
% %                 file_name = sprintf('%sDP_par_simulation%d_choice%d_capacity%d_nopurch%d_diffa1%d_noExtra%d.mat',folder,nsim,choice,c,a0,a1,noExtra);
% %                 save(file_name,'V','Vroot','t_training','policy2','a0','a1','a2','a3','b1','b2','b3','tau1','tau2','tau3','c','T','rev_nt','rev','policy','tickets_n','seats2_n','seats3_n');
% %                 filename = 'Results_DP.xlsx';
% %                 writematrix('',filename,'Sheet',sheet,'Range','A'+string(i));
% %                 writematrix(c,filename,'Sheet',sheet,'Range','B'+string(i));
% %                 writematrix(Vroot,filename,'Sheet',sheet,'Range','C'+string(i));
% %                 writematrix(t_training,filename,'Sheet',sheet,'Range','D'+string(i));
% %                 writematrix(rev,filename,'Sheet',sheet,'Range','E'+string(i));
% %                 writematrix(t_sim,filename,'Sheet',sheet,'Range','F'+string(i));
% %                 writematrix(nsim,filename,'Sheet',sheet,'Range','G'+string(i));
% %                 writematrix(a0,filename,'Sheet',sheet,'Range','H'+string(i));
% % %                 writematrix(tau1,filename,'Sheet',sheet,'Range','I'+string(i));
% %                 writematrix(stderr,filename,'Sheet',sheet,'Range','J'+string(i));
% %                 writematrix(mean(tickets_n),filename,'Sheet',sheet,'Range','K'+string(i)+':N'+string(i));
% %                 writematrix(mean(seats2_n),filename,'Sheet',sheet+1,'Range',string(i)+':'+string(i));
% %                 writematrix(mean(seats3_n),filename,'Sheet',sheet+2,'Range',string(i)+':'+string(i));
% %                 i = i+1;
% %             end
% %         end
% %     end
% % end
% 
% % %% PriceSensitivity %%
% % for noExtra = [1] %=0: selling product 3, =1: not selling product 3 = the price of product 3 = inf
% %     for choice = [0,1,2]
% %         if noExtra == 0
% %             i = 14;
% %         else
% %             i = 25;
% %         end
% %         if choice == 0
% %             sheet = 1;
% %         elseif choice == 1
% %             sheet = 4;
% %         elseif choice == 2
% %             sheet = 7;
% %         end
% %         for c = [8] % Never chose an odd number!
% %             T = 2*c;
% % 
% %             % Generate values in MATLAB
% %             if choice == 0
% %                 [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_homo_seats(c);
% %                 folder = 'SensitivityAnalysisResults/LevelDissim23/HOMOG/';
% %             elseif choice == 1
% %                 [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_aisle_window_seats(c);
% %                 folder = 'SensitivityAnalysisResults/LevelDissim23/AW/';
% %             elseif choice == 2
% %                 [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_increasing_seats(c);
% %                 folder = 'SensitivityAnalysisResults/LevelDissim23/HET/';
% %             end
% %             disp([a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3]);
% %             bmat = [0.5*[b1;b2;b3],0.6*[b1;b2;b3],0.7*[b1;b2;b3],0.8*[b1;b2;b3],0.9*[b1;b2;b3],[b1;b2;b3],1.1*[b1;b2;b3],1.2*[b1;b2;b3],1.3*[b1;b2;b3],1.4*[b1;b2;b3],1.5*[b1;b2;b3]];
% % 
% %             for b = bmat
% %                 b1 = b(1);
% %                 b2 = b(2);
% %                 b3 = b(3);
% % 
% %                 % Initialize V table with zeros: V[seat arrangement,reservation,time]
% %                 V = zeros(2^c + 1,c + 2,T + 1);
% %                 V(:,1,:) = -9999;
% %                 V(1,:,:) = -9999;
% %                 Vhelp = zeros(2^c + 1,c + 2);
% %                 Vhelp(:,1) = -9999;
% %                 Vhelp(1,:) = -9999;
% %                 policy = zeros(2^c + 1,c + 2,1+2*c,T+1);
% %                 policyhelp = zeros(2^c + 1,c+2,1+2*c);
% % 
% %                 % Learning process
% %                 Y = sum(c)+2;
% %                 tstart = tic;
% %                 for t = 2:T+1
% %                     disp(t-1);
% %                     parfor D = 3:(2^c+1)
% %                         x = DtoB0(D,c);
% %                         for y = 3:Y 
% %                             if y>sum(x)+2
% %                                 Vhelp(D,y) = -9999;
% %                             else
% %                                 [ofv,r1,r2j,r3j] = opt_DP(V,c,x,y,t,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,noExtra);
% %                                 Vhelp(D,y) = ofv;
% %                                 policyhelp(D,y,:) = [r1;r2j;r3j];
% %                             end
% %                         end
% %                     end
% %                     policy(:,:,:,t) = policyhelp;
% %                     V(:,:,t) = Vhelp;
% %                 end
% %                 t_training = toc(tstart);
% % 
% %                 x = ones(c,1);
% %                 y = sum(x) + 2;
% %                 t = T+1;
% %                 disp(V(BtoD0(x),y,t));
% %                 Vroot = V(BtoD0(x),y,t);
% %                 policy2 = policy;
% % 
% %                 % Get random decision variables
% %                 rng(987654);
% %                 nsim = 100;
% %                 dec = sample_dec(T+1,nsim);
% %                 tstart = tic;
% %                 [rev_nt,policy,state,tickets_n,seats2_n,seats3_n] = simulation(nsim,dec,T,c,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,V,noExtra);
% %                 t_sim = toc(tstart);
% %                 stderr = std(sum(rev_nt,2));
% %                 rev = mean(sum(rev_nt,2))
% % 
% %                 %save data
% %                 file_name = sprintf('%sDP_par_simulation%d_choice%d_capacity%d_nopurch%d_pricesenb1%_noExtra%d.mat',folder,nsim,choice,c,a0,b1,noExtra);
% %                 save(file_name,'V','Vroot','t_training','policy2','a0','a1','a2','a3','b1','b2','b3','tau1','tau2','tau3','c','T','rev_nt','rev','policy','tickets_n','seats2_n','seats3_n');
% %                 filename = 'Results_DP.xlsx';
% %                 writematrix('',filename,'Sheet',sheet,'Range','A'+string(i));
% %                 writematrix(c,filename,'Sheet',sheet,'Range','B'+string(i));
% %                 writematrix(Vroot,filename,'Sheet',sheet,'Range','C'+string(i));
% %                 writematrix(t_training,filename,'Sheet',sheet,'Range','D'+string(i));
% %                 writematrix(rev,filename,'Sheet',sheet,'Range','E'+string(i));
% %                 writematrix(t_sim,filename,'Sheet',sheet,'Range','F'+string(i));
% %                 writematrix(nsim,filename,'Sheet',sheet,'Range','G'+string(i));
% %                 writematrix(a0,filename,'Sheet',sheet,'Range','H'+string(i));
% % %                 writematrix(tau1,filename,'Sheet',sheet,'Range','I'+string(i));
% %                 writematrix(stderr,filename,'Sheet',sheet,'Range','J'+string(i));
% %                 writematrix(mean(tickets_n),filename,'Sheet',sheet,'Range','K'+string(i)+':N'+string(i));
% %                 writematrix(mean(seats2_n),filename,'Sheet',sheet+1,'Range',string(i)+':'+string(i));
% %                 writematrix(mean(seats3_n),filename,'Sheet',sheet+2,'Range',string(i)+':'+string(i));
% %                 i = i+1;
% %             end
% %         end
% %     end
% % end
% 
