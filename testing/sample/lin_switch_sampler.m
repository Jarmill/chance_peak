function [out_sim] = lin_switch_sampler(x0, Tmax, dt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
            if nargin < 5
                curr_event = @obj.loc.supp_event;                
            end
            
            %copy over from peak/sample_cont.m
            %main solving loop
            x0_curr = x0;
            Nperiod = ceil(Tmax/dt);
            % time_accum = zeros(1, Nperiod);
            time_accum = 0;
            % x_accum = zeros(2, Nperiod);
            x_accum = reshape(x0, 1, []);
            time_index= 0;

            time_total = 0;
            
            system_choice = [];
            time_breaks = 0;
            k = 1;

            FINE = 1;

    

%             Antithetic = 1;
            NTrials = 1;
            dt = 1e-3;

            mu = 0.5;

            curr_sys_ind = (rand(1)>0.5);

%             simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', NTrials,...
%       'Antithetic', Antithetic);

            f1 = @(t, x) [-5, -4; -1, -2]*x/2;
            g1 = @(t, x) [0; 0.5*x(2)]/2;


            f2 = @(t, x) [-2, -4; 5, -2]*x/2;
            g2 = @(t, x) g1(t, x);
            
            while time_total < Tmax
    
                %choose a possible switching system that is admissible for current time/state               

                %for how long should this uncertain system be sampled?

                time_track = max(exprnd(mu, 1, 1), dt);



                %do not exceed Tmax
                time_track_trunc = min(time_track, Tmax - time_total);

                if (Tmax - time_total) <= dt
                    break
                end
%                 sde_curr
                %simulate the current system
                if FINE
                    Nperiod_curr = ceil(time_track_trunc/dt);
                    curr_ode_options =   odeset('Events',curr_event, 'RelTol', 1e-7, ...
                                                  'AbsTol', 1e-8, 'MaxStep', 0.01);
                else
                    curr_ode_options =   odeset('Events',curr_event);
                end


                %MATLAB Financial toolbox for sampling
%                 sde_curr = sde(curr_f, curr_g, 'StartState', x0_curr);    % dX = F(t,X)dt + G(t,X)dW
%                 [x_curr,time_curr] = simByEuler(sde_curr, Nperiod, 'DeltaTime', dt, 'NTrials', NTrials);
                
                %SDETOOLS for sampling

                if curr_sys_ind == 0
                    curr_f = f1;
                    curr_g = g1;
                else
                    curr_f = f2;
                    curr_g = g2;
                end

                curr_event = @(t, x) box_event(t, x, 3);
                %  
                curr_sde_options = sdeset('EventsFun', curr_event,'SDEType', 'Ito', ...
                    'DiagonalNoise', 'yes', 'ConstGFUN', 'yes');
                [x_curr,W,TE,YE,WE,IE] = sde_euler(curr_f,curr_g,0:dt:time_track_trunc,...
                    x0_curr,curr_sde_options);
                
                time_curr = (1:(size(x_curr, 1)-1))*dt;
                time_end = time_curr(end);
                %save current trajectory
                %check indices/dimensions
                x0_curr = x_curr(end, :)';
                x_accum = [x_accum; x_curr(2:end, :)];
                time_accum = [time_accum, time_curr + time_total];

                time_total = time_total + time_curr(end);

                %time dependent uncertainty
                system_choice = [system_choice; curr_sys_ind];  %system switching
                time_breaks = [time_breaks; time_total];

                time_index= [time_index; k*ones(length(time_curr), 1)];
                k = k + 1;
                curr_sys_ind = 1-curr_sys_ind;
            end
            
%             %simulate the trajectory
%             curr_ode_options = odeset('Events',curr_event, 'RelTol', 1e-7, ...
%                                       'AbsTol', 1e-8, 'MaxStep', 0.01);
%             
            %package up the output            
            out_sim = struct;

            %trajectories
            out_sim.t = time_accum;
            out_sim.x = x_accum;
            out_sim.system = system_choice;
            out_sim.time_breaks = time_breaks;
           
end

function [event_eval, terminal, direction] = box_event(t, x, box)
%BOX_EVENT Event function for switch_sim used when simulating systems
%   return the event function eventFcn. Yes this is tricky
%is always true.

    %stop integrating when the system falls outside box
    %
    event_eval = all(abs(x) <= box);
    terminal = 1;
    direction = 0;
end
