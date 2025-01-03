for noExtra = [0] %=0: selling product 3, =1: not selling product 3 = the price of product 3 = inf
    for choice = [0]
        if noExtra == 0
            i = 28;
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
            a3val = 0:0.1:5;
            a3mat = repmat(a3val,c,1);
            for a3 = a3mat
                a3 = a3';
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
                rng(9876543);
                nsim = 100;
                dec = sample_dec(T+1,nsim);
                tstart = tic;
                [rev_nt,policy,state,tickets_n,seats2_n,seats3_n] = simulation(nsim,dec,T,c,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,V,noExtra);
                t_sim = toc(tstart);
                stderr = std(sum(rev_nt,2));
                rev = mean(sum(rev_nt,2))

                %save data
                file_name = sprintf('%sDP_par_simulation%d_choice%d_capacity%d_noExtra%d_a3_%d.mat',folder,nsim,choice,c,noExtra,a3(1));
                save(file_name,'dec','V','Vroot','t_training','SLF','policy2','a0','a1','a2','a3','b1','b2','b3','tau1','tau2','tau3','c','T','rev_nt','rev','policy','tickets_n','seats2_n','seats3_n');
                filename = 'Results_DP.xlsx';
                sheet = 15; 
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
%                 writematrix(mean(seats2_n),filename,'Sheet',sheet+1,'Range',string(i)+':'+string(i));
%                 writematrix(mean(seats3_n),filename,'Sheet',sheet+2,'Range',string(i)+':'+string(i));
                i = i+1;
            end
        end
    end
end
