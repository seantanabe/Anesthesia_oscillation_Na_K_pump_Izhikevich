%% U_CC, stimuli Wake J 0.005 vs. SO J 0.0005

clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
addpath('F:\20190624_Wisconsin\violin_plot')

OUT_states = cell(2,4);
for state = 1:2
    
    if state == 1
%         Analysis  = load('20201125T184545_Izhi_run_J0005_stimuli.mat') % after stimuli 100s J 0.005 (Wake)
%         Analysis  = load('20201223T211706_Izhi_run_J0005_stimuli_v2.mat') % after stimuli 100s J 0.005 (Wake)
%         Analysis  = load('20201228T235903_Izhi_run_J0005_stimuli_randi_v3.mat') % after rand-stimuli 100s J 0.005 (Wake)
        Analysis  = load('20201231T192605_Izhi_run_J0005_stimuli_randi_v4.mat') % after rand-stimuli 100s J 0.005 (Wake)
    elseif state == 2
%         Analysis  = load('20201125T183116_Izhi_run_J00005_stimuli.mat') % after stimuli 100s J 0.0005 (SO)
%         Analysis  = load('20201224T005147_Izhi_run_J00005_stimuli_v2.mat') % after stimuli 100s J 0.0005 (SO)
%         Analysis  = load('20201229T031229_Izhi_run_J00005_stimuli_randi_v3.mat') % after rand-stimuli 100s J 0.0005 (SO)
        Analysis  = load('20201231T225045_Izhi_run_J00005_stimuli_randi_v4.mat') % after rand-stimuli 100s J 0.0005 (SO)
    end
    
    h      = Analysis.param.h;
    N      = Analysis.param.N;
    trials = Analysis.initial.T_sec ;
    
    stim_onset    = 250;
    segment       = 100; %200; %ms
    segmentPre    = 200;
    convFilter    = gausswin(4/h); %%%%%%%%%%%%%%%%%%%%%%
    
    firings_store = Analysis.trend.firings_store;
    firings_store = squeeze(sum(firings_store,1));
    firings_store = firings_store(:,1:(1000/(2*h))) + firings_store(:,((1000/(2*h))+1):end);
    firings_store = conv2(firings_store,convFilter','same');
    
    firings_store_prestim  = firings_store(:,((stim_onset-segmentPre-1)/h):((stim_onset-1)/h));
    firings_store_poststim = firings_store(:,((stim_onset+1)/h):((stim_onset+segment+1)/h));
    
    
    
    conv_sig = (firings_store_poststim' > prctile(firings_store_prestim',95))'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARMA
    conv_sig = [zeros(N,excl_len) conv_sig];
    conv_sig = [conv_sig zeros(N,excl_len)];
    for i = 1:(size(conv_sig,2)-excl_len)
        for n = 1:size(conv_sig,1)
            if sum(conv_sig(n,[i (i+excl_len)]),2) == 0
                conv_sig(n,i:(i+excl_len)) = 0;
            end
        end
    end
    conv_sig = conv_sig(:,(excl_len+1):end);
    conv_sig = conv_sig(:,1:(end-excl_len-1));
    
    conv_sig_firing = conv_sig.*firings_store_poststim(:,1:end-1);
    
    cenmass = zeros(N,1);
    for i = 1:N
        cenmass(i) = sum(conv_sig_firing(i,:).*(h:h:(segment)))/sum(conv_sig_firing(i,:));
    end
    %     cenmass = zeros(N,1);
    %     for i = 1:N
    %         cenmass(i) = sum(firings_store_poststim(i,:).*(h:h:(segment+h)))/sum(firings_store_poststim(i,:));
    %     end
    %
    excl_n = sum(conv_sig,2) == 0 ;
    sig_neuron = sum(conv_sig,2) ~= 0 ;
    
    OUT_states{state,1} = firings_store_poststim';
    OUT_states{state,2} = conv_sig;
    OUT_states{state,3} = excl_n;
    OUT_states{state,4} = cenmass;
    
    OUT_states{state,5} = [ones(800,1); zeros(200,1)];
    OUT_states{state,6} = Analysis.param.D;
    %
    if state == 1 % after stimuli 100s J 0.005 (Wake)
        cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
        save(['20201231_Izhi_run_J0005_W_stimuli_cenmass' ],'cenmass');
        save(['20201231_Izhi_run_J0005_W_stimuli_sig_neuron' ],'sig_neuron');
        save(['20201231_Izhi_run_J0005_W_stimuli_firings_store_prestim' ],'firings_store_prestim');
        save(['20201231_Izhi_run_J0005_W_stimuli_firings_store_poststim' ],'firings_store_poststim');
    elseif state == 2 % after stimuli 100s J 0.0005 (SO)
        cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
        save(['20201231_Izhi_run_J00005_SO_stimuli_cenmass' ],'cenmass');
        save(['20201231_Izhi_run_J00005_SO_stimuli_sig_neuron' ],'sig_neuron');
        save(['20201231_Izhi_run_J00005_SO_stimuli_firings_store_prestim' ],'firings_store_prestim');
        save(['20201231_Izhi_run_J00005_SO_stimuli_firings_store_poststim' ],'firings_store_poststim');
        
    end
    
end