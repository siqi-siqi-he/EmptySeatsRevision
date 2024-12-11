for noExtra = [1] %=0: selling product 3, =1: not selling product 3 = the price of product 3 = inf
    for choice = [0,1,2]
        if noExtra == 0
            i = 2;
        else
            i = 10;
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
                V(:,1,:) = -9999;
                V(1,:,:) = -9999;
                Vhelp = zeros(2^c + 1,c + 2);
                Vhelp(:,1) = -9999;
                Vhelp(1,:) = -9999;

                % Learning process
                Y = sum(c)+2;
                tstart = tic;
                for t = 2:T+1
                    disp(t-1);
                    for D = 3:(2^c+1)
                        x = DtoB0(D,c);
                        for y = 3:Y %Y(D) does not work due to parfor, thus another if clause
                            if y>sum(x)+2
                                Vhelp(D,y) = -9999;
                            else
                                % Optimize prices on expectation value
                                ofv = opt_DP(V,c,x,y,t,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,noExtra);
                                Vhelp(D,y) = ofv;
                            end
                        end
                    end
                    V(:,:,t) = Vhelp;
                end
                t_training = toc(tstart);

                x = ones(c,1);
                y = sum(x) + 2;
                t = T+1;
                disp(V(BtoD0(x),y,t));
                Vroot = V(BtoD0(x),y,t);

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
                save(file_name,'V','Vroot','t_training','a0','a1','a2','a3','b1','b2','b3','tau1','tau2','tau3','c','T','rev_nt','rev','policy','tickets_n','seats2_n','seats3_n');
                filename = 'Results_DP.xlsx';
                writematrix('',filename,'Sheet',sheet,'Range','A'+string(i));
                writematrix(c,filename,'Sheet',sheet,'Range','B'+string(i));
                writematrix(Vroot,filename,'Sheet',sheet,'Range','C'+string(i));
                writematrix(t_training,filename,'Sheet',sheet,'Range','D'+string(i));
                writematrix(rev,filename,'Sheet',sheet,'Range','E'+string(i));
                writematrix(t_sim,filename,'Sheet',sheet,'Range','F'+string(i));
                writematrix(nsim,filename,'Sheet',sheet,'Range','G'+string(i));
                writematrix(a0,filename,'Sheet',sheet,'Range','H'+string(i));
                writematrix(0,filename,'Sheet',sheet,'Range','I'+string(i));
                writematrix(stderr,filename,'Sheet',sheet,'Range','J'+string(i));
                writematrix(mean(tickets_n),filename,'Sheet',sheet,'Range','K'+string(i)+':N'+string(i));
                writematrix(mean(seats2_n),filename,'Sheet',sheet+1,'Range',string(i)+':'+string(i));
                writematrix(mean(seats3_n),filename,'Sheet',sheet+2,'Range',string(i)+':'+string(i));
                i = i+1;
            end
        end
    end
end
