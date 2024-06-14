
clear all
rand('seed',1);
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
addpath('E:\poissonTutorial')

initial = [];
param   = [];

param.N       = 1000;
param.density = 100; % N/mm2
param.N_inh   = round(param.N*0.2);
param.N_ex    = param.N - param.N_inh;

% letter stimuli
[initial.connect, initial.N_coord]    = connect_plate_periodic(param.N,param.density); % row to, col from
initial.stim_coord                    = [1.5; 1.5];
initial.sigma_stim                    = 0.2;
[initial.letter,initial.stim_neurons] = letter_gauss(param,initial);
param.stimuliHz                       = 0;

initial.v         = -70*ones(param.N,1);
initial.u         = 0.2*(initial.v);
initial.g_max     = 0.045;
% initial.g_max     = 0.1;
% initial.g_peak    = [initial.connect(:,1:param.N_ex).*(initial.g_max) ...
%                      initial.connect(:,(param.N_ex+1):param.N)*1]; %05];
% initial.g_peak    = [initial.connect(:,1:param.N_ex).*(initial.g_max) ...
%                      initial.connect(:,(param.N_ex+1):param.N)*1]; %05];
% initial.g_peak    = [initial.connect(:,1:param.N_ex).*(initial.g_max*rand(param.N,param.N_ex)) ...
%                      initial.connect(:,(param.N_ex+1):param.N)*0.05];
initial.g_peak   = [initial.connect(:,1:param.N_ex).*(initial.g_max/10) ...
    initial.connect(:,(param.N_ex+1):param.N)*0.05];
% initial.g_peak   = [initial.connect(:,1:param.N_ex).*(initial.g_max/10) ...
%                     initial.connect(:,(param.N_ex+1):param.N)*0.015];
initial.g         = zeros(param.N,param.N);
initial.T_sec     = 10;
initial.C_ATP     = ones(param.N,1)*2 ;
initial.C_Na      = ones(param.N,1)*25 ;
initial.z         = zeros(param.N,1);
initial.g_Na_peak = 50*(0.1);
initial.g_Na      = zeros(param.N,1);

param.stdp1        = 1;
param.display      = 1;
param.h            = 0.1 ;
param.A_pos        = 0.005;
param.noiseHz      = 0; %%%.5;
param.noise_amp_ex = 0; %%%%%%%%30;
param.noise_amp_in = 0 ;
param.stimuli_onset_1s   = zeros(1,1000/(param.h));
param.J_ATP        = 0.1 ; %0.1;  %%%%%%%%%%%%%%%%%
param.I_app_DC     = [ones(param.N_ex,1).*1.8*22; ... %17; ...
    ones(param.N_inh,1).*0.6*22]; %17];
param.re           =rand(param.N_ex,1);
param.ri           =rand(param.N_inh,1);
param.A            =[0.02*ones(param.N_ex,1); 0.02+0.08*(param.ri)];
param.B            =[0.2*ones(param.N_ex,1);  0.25-0.05*(param.ri)];
param.C            =[-65+15*(param.re).^2;    -65*ones(param.N_inh,1)];
param.D            =[8-6*(param.re).^2;       2*ones(param.N_inh,1)];

[initial, trend] = Izhi_ATP_run_v6(initial,param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_Izhi_run' ],'-v7.3');


%% find groups

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005
clear trend

param.anchor_width   = 3;
param.T_group        = 50; %ms
param.min_group_path = 2;
param.display        = 0;

[groups] = group_find_v1(initial,param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_groups_v2' ],'-v7.3');
% 20201224T104121_groups_v2

%% group analysis plot; size, path, example

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat') %

group_n = size(groups,2);
size_dt = zeros(1,group_n);
path_dt = zeros(1,group_n);
for i = 1:size(groups,2)
    size_dt(1,i) = groups{i}.size;
    path_dt(1,i) = groups{i}.longest_path;
end

figure
histogram(size_dt,'FaceColor', 'k')
xlabel('group size'); ylabel('count')

figure
histogram(path_dt,'FaceColor', 'k')
xlabel('longest path'); ylabel('count')

path_4_ind = find(path_dt == 4);
ex_gr = 639;
gr = groups{ex_gr};
h = param.h;
T_group = param.T_group;
N = param.N;

figure
[fire_I,fire_J] = find(gr.firings);
fire_J = fire_J*h;
scatter(fire_J,fire_I, 50,'ok');
xlim([0 25]) %T_group])
ylim([1 N])
hold on
for i_conn = 1:size(gr.gr_conn,1)
    plot(gr.gr_conn(i_conn,[2 3 5]),gr.gr_conn(i_conn,[1 4 4]), 'k')
end
hold off
ylabel('neurons')
xlabel('ms')
title(['group size ' num2str(gr.size) ', longest path ' num2str(gr.longest_path)])
%% group analysis; plot group in space

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
load('20201224T104121_groups_v2.mat'); %
% Data = load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat'); % after 100s J 0.005
Data = load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat'); % after 100s J 0.0005


% Analysis = load('20210118T020811_group_active_find_J0005.mat') % W
Analysis = load('20210118T021757_group_active_find_J00005.mat') % SO

% firings_store    = Data.trend.firings_store;
firings_store    = Analysis.Data.trend.firings_store;
N                = param.N;
h                = param.h;
N_coord = initial.N_coord;

n_ind_s = round(1000/h);
n_trial = size(firings_store,1);

firing_long = zeros(N, n_ind_s*n_trial);
for s_i = 1:n_trial
    firing_long(:,((s_i-1)*n_ind_s+1):(s_i*n_ind_s)) = squeeze(firings_store(s_i,:,:));
end


gr_size_all = zeros(1,size(groups,2));
for gr_i = 1:size(groups,2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['find group ' num2str(gr_i) ' active'])
    gr_size_all(gr_i)             = groups{gr_i}.longest_path; %%%%%%%%%not size
end
[M,I] = max(gr_size_all) 
find(gr_size_all == M)

N                = param.N;
h                = param.h;
gr_active_jit    = 1;

jit_kernel = ones(1,round(gr_active_jit/h)*2);
jit_ind    = round(gr_active_jit/h)*2;
group_n    = size(groups,2);

gr_i = 1763 %%%%%%%%%%%%%%%%%%% gr_1
gr_kernel_orig      = groups{gr_i}.firings;
gr_anchor_neurons   = groups{gr_i}.anchor_neurons;
[gr_ind_spikes_I,gr_ind_spikes_J]  = find(gr_kernel_orig);
for i = 1:length(gr_anchor_neurons)
    tmp_ind = find((gr_ind_spikes_I == gr_anchor_neurons(i)));
    gr_kernel_orig(gr_ind_spikes_I(tmp_ind(1)),gr_ind_spikes_J(tmp_ind(1))) = 0;
end
gr_kernel   = conv2(gr_kernel_orig,jit_kernel,'full');
ind_ker = find(sum(gr_kernel,2) > 0);
gr_1 = ind_ker; %%%%%%%%%%%%% gr_1

gr_i = 1680   %%%%%%%%%%%%%%%%%% gr_2
gr_kernel_orig      = groups{gr_i}.firings;
gr_anchor_neurons   = groups{gr_i}.anchor_neurons;
[gr_ind_spikes_I,gr_ind_spikes_J]  = find(gr_kernel_orig);
for i = 1:length(gr_anchor_neurons)
    tmp_ind = find((gr_ind_spikes_I == gr_anchor_neurons(i)));
    gr_kernel_orig(gr_ind_spikes_I(tmp_ind(1)),gr_ind_spikes_J(tmp_ind(1))) = 0;
end
gr_kernel   = conv2(gr_kernel_orig,jit_kernel,'full');
ind_ker = find(sum(gr_kernel,2) > 0);
gr_2 = ind_ker; %%%%%%%%%%%%%%%% gr_2

gr_i = 2500   %%%%%%%%%%%%%%%%%% gr_3
gr_kernel_orig      = groups{gr_i}.firings;
gr_anchor_neurons   = groups{gr_i}.anchor_neurons;
[gr_ind_spikes_I,gr_ind_spikes_J]  = find(gr_kernel_orig);
for i = 1:length(gr_anchor_neurons)
    tmp_ind = find((gr_ind_spikes_I == gr_anchor_neurons(i)));
    gr_kernel_orig(gr_ind_spikes_I(tmp_ind(1)),gr_ind_spikes_J(tmp_ind(1))) = 0;
end
gr_kernel   = conv2(gr_kernel_orig,jit_kernel,'full');
ind_ker = find(sum(gr_kernel,2) > 0);
gr_3 = ind_ker; %%%%%%%%%%%%%%%% gr_3

figure
scatter(N_coord(1,gr_1),N_coord(2,gr_1),100,'xr')
hold on
scatter(N_coord(1,gr_2),N_coord(2,gr_2),100,'ob')
hold on
scatter(N_coord(1,gr_3),N_coord(2,gr_3),100,'*g')
hold on
scatter(N_coord(1,:),N_coord(2,:),10,'.k')
hold on
legend('example group 1','example group 2','example group 3')
xlabel('mm')
ylabel('mm')


% x_length = [0 1]*(1000/h)+1
% x_length = [0.32 0.335]*(1000/h)+1 % gr_1
% x_length = [0.345 0.36]*(1000/h)+1 % gr_3
x_length = [0.375 0.39]*(1000/h)+1   % gr_2

firing_set = sum(firing_long(:,x_length(1):x_length(2)),2);

max(sum(firing_long,2))

figure
scatter(N_coord(1,:),N_coord(2,:),100,firing_set,'.')
colormap(flipud(hot))
caxis([0 3])
colorbar

%% group analysis; activation, Wake (J0.005) 

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat'); %
% Data = load('20201031T224652_Izhi_run_J01.mat') %%%%%%%%%%%%%%%%
Data = load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat'); % after 100s J 0.005
% load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works

Data.param.gr_active_jit   = 1;

trend = Data.trend;
param = Data.param;
[groups] = group_active_find_v4(groups,Data.trend,Data.param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J0005' ],'groups','Data','-v7.3');

%% group analysis; activation, SO (J0.0005)

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat'); %
% load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
Data = load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat'); % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works

Data.param.gr_active_jit   = 1;

trend = Data.trend;
param = Data.param;
[groups] = group_active_find_v4(groups,Data.trend,Data.param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J00005' ],'groups','Data','-v7.3');

%% group analysis (CONTROL jitter); activation, Wake (J0.005) 

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat'); %
% Data = load('20201031T224652_Izhi_run_J01.mat') %%%%%%%%%%%%%%%%
Data = load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat'); % after 100s J 0.005
% load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works
% firings_store_orig = Data.trend.firings_store;
%%%%%%% jitter
trend                 = Data.trend;
Data.param.jit_firing = 2*20; %2*0.2;
param                 = Data.param;

Data.trend.firings_store = firing_jitter_v1(trend, param);
% firings_store_new = Data.trend.firings_store;
%%%%%%

% sum(sum(sum(abs(firings_store_new - firings_store_orig))))

Data.param.gr_active_jit   = 1;

trend = Data.trend;
param = Data.param;
[groups] = group_active_find_v4(groups,Data.trend,Data.param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J0005_jitter' ],'groups','Data','-v7.3');

%% group analysis (CONTROL jitter); activation, SO (J0.0005)

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat'); %
% load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
Data = load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat'); % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works

%%%%%%% jitter
trend                 = Data.trend;
Data.param.jit_firing = 2*40; %2*0.2;
param                 = Data.param;

Data.trend.firings_store = firing_jitter_v1(trend, param);

%%%%%%

Data.param.gr_active_jit   = 1;

trend = Data.trend;
param = Data.param;
[groups] = group_active_find_v4(groups,Data.trend,Data.param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J00005_jitter' ],'groups','Data','-v7.3');

%% group analysis (CONTROL jitter); activation, Wake (J0.005) x 20 rand

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat'); %
% Data = load('20201031T224652_Izhi_run_J01.mat') %%%%%%%%%%%%%%%%
Data = load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat'); % after 100s J 0.005
% load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works

trend_orig                 = Data.trend;
param_orig                 = Data.param;
groups_orig                = groups;

param_orig.jit_firing = 2*2;
Data.param.jit_firing = 2*2; %2*0.2;

groups_perm = cell(1,length(groups));
for r = 1:20
    Data.trend.firings_store = firing_jitter_v1(trend_orig, param_orig);
    
    Data.param.gr_active_jit   = 1;
    
    [groups] = group_active_find_v4(groups_orig,Data.trend,Data.param);
    
    for i = 1:length(groups)
        groups_perm{i} = [groups_perm{i} groups{i}.active_per_s];
    end
end

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J0005_jitter_2ms_groups_perm' ],'groups_perm','-v7.3');

%% group analysis (CONTROL jitter); activation, SO (J0.0005) x 20 rand

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat'); %
% load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
Data = load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat'); % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works

trend_orig                 = Data.trend;
param_orig                 = Data.param;
groups_orig                = groups;

param_orig.jit_firing = 2*2;
Data.param.jit_firing = 2*2; %2*0.2;

groups_perm = cell(1,length(groups));
for r = 1:20
    Data.trend.firings_store = firing_jitter_v1(trend_orig, param_orig);
    
    Data.param.gr_active_jit   = 1;
    
    [groups] = group_active_find_v4(groups_orig,Data.trend,Data.param);
    
    for i = 1:length(groups)
        groups_perm{i} = [groups_perm{i} groups{i}.active_per_s];
    end
end

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J00005_jitter_2ms_groups_perm' ],'groups_perm','-v7.3');

%% group analysis (CONTROL FlipRast) ; activation, Wake (J0.005)

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat') %
% Data = load('20201031T224652_Izhi_run_J01.mat') %%%%%%%%%%%%%%%%
Data = load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
% load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works

Data.param.gr_active_jit   = 1;

Data.trend.firings_store = flip(Data.trend.firings_store,3);
Data.trend.firings_store = flip(Data.trend.firings_store,1);

% [groups] = group_active_find_v3(groups,Data.trend,Data.param);
% [groups] = group_active_find_v2(groups,Data.trend,Data.param);
[groups] = group_active_find_v4(groups,Data.trend,Data.param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J0005_FlipRast' ],'groups','Data','-v7.3');

%% group analysis (CONTROL FlipRast); activation, SO (J0.0005)

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat') %
% load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
Data = load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works

Data.param.gr_active_jit   = 1;

Data.trend.firings_store = flip(Data.trend.firings_store,3);
Data.trend.firings_store = flip(Data.trend.firings_store,1);

% [groups] = group_active_find_v3(groups,Data.trend,Data.param);
% [groups] = group_active_find_v2(groups,Data.trend,Data.param);
[groups] = group_active_find_v4(groups,Data.trend,Data.param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J00005_FlipRast' ],'groups','Data','-v7.3');


%% group analysis (CONTROL FlipRast ShuffN) ; activation, Wake (J0.005)

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat') %
% Data = load('20201031T224652_Izhi_run_J01.mat') %%%%%%%%%%%%%%%%
Data = load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
% load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works

Data.param.gr_active_jit   = 1;

Data.trend.firings_store = flip(Data.trend.firings_store,3);
Data.trend.firings_store = flip(Data.trend.firings_store,1);

Data.trend.firings_store = Data.trend.firings_store(:,randperm(Data.param.N),:);

% [groups] = group_active_find_v3(groups,Data.trend,Data.param);
% [groups] = group_active_find_v2(groups,Data.trend,Data.param);
[groups] = group_active_find_v4(groups,Data.trend,Data.param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J0005_FlipRast_ShuffN' ],'groups','Data','-v7.3');

%% group analysis (CONTROL FlipRast ShuffN); activation, SO (J0.0005)

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat') %
% load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
Data = load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005

Data.trend.firings_store = Data.trend.firings_store((end-60+1):end,:,:); %% memory short, check later, but this works

Data.param.gr_active_jit   =1;

Data.trend.firings_store = flip(Data.trend.firings_store,3);
Data.trend.firings_store = flip(Data.trend.firings_store,1);

Data.trend.firings_store = Data.trend.firings_store(:,randperm(Data.param.N),:);

% [groups] = group_active_find_v3(groups,Data.trend,Data.param);
% [groups] = group_active_find_v2(groups,Data.trend,Data.param);
[groups] = group_active_find_v4(groups,Data.trend,Data.param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_find_J00005_FlipRast_ShuffN' ],'groups','Data','-v7.3');

%% group analysis; activation plot W and SO

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
addpath('F:\20190624_Wisconsin\violin_plot')

normalization = 0; % 0: no norm, 1: norm FR all, 2: norm FR gr

% wake
Analysis = load('20210126T230008_group_active_find_J0005.mat') % W v4 filter2 anchor exclude 50%
% Perm = load('20210126T115509_group_active_find_J0005_jitter_5ms_groups_perm.mat') % W v4 filter2 anchor exclude 50% perm jitter 5ms
Perm = load('20210127T034107_group_active_find_J0005_jitter_1ms_groups_perm.mat') % W v4 filter2 anchor exclude 50% perm jitter 1ms

firings_store   = Analysis.Data.trend.firings_store ;
n_trial         = size(firings_store,1);
FR = sum(sum(firings_store,1),3)/n_trial;

group_n   = size(Analysis.groups,2);
active_dt             = zeros(1,group_n);
FR_gr  = zeros(1,group_n);
CI95 = zeros(1,group_n);
for i = 1:group_n
    active_dt(1,i) = Analysis.groups{i}.active_per_s;
    tmp_neurons = Analysis.groups{i}.neurons;
    FR_gr(i) = mean(sum(sum(firings_store(:,tmp_neurons,:),1),3)/n_trial);
    CI95(i) = prctile(Perm.groups_perm{i},0.95);
end
active_dt_norm_FR_all  = active_dt/sum(FR);
active_dt_norm_FR_gr = active_dt./FR_gr;

CI95_norm_FR_all  = CI95/sum(FR);
CI95_norm_FR_gr = CI95./FR_gr;

if normalization == 0
    active_W_dt     = active_dt;
    CI95_W = (active_W_dt >= CI95);
elseif normalization == 1
    active_W_dt     = active_dt_norm_FR_all;
    CI95_W = (active_W_dt >= CI95_norm_FR_all);
elseif normalization == 2
    active_W_dt     = active_dt_norm_FR_gr;
    CI95_W = (active_W_dt >= CI95_norm_FR_gr);
end
sum(CI95_W)/group_n 

figure
histogram(CI95)
hold on
histogram(active_W_dt)

figure
scatter(CI95, active_W_dt,'.k')
hold on
plot([0 1400], [0 1400])

% SO
Analysis   = load('20210126T230848_group_active_find_J00005.mat') % SO v4 filter2 exc anchor %50
% Perm = load('20210126T140716_group_active_find_J00005_jitter_5ms_groups_perm.mat') % SO v4 filter2 anchor exclude 50% perm jitter 5ms
Perm = load('20210127T053920_group_active_find_J00005_jitter_1ms_groups_perm.mat') % SO v4 filter2 anchor exclude 50% perm jitter 1ms

firings_store   = Analysis.Data.trend.firings_store ;
n_trial         = size(firings_store,1);
FR = sum(sum(firings_store,1),3)/n_trial;

group_n   = size(Analysis.groups,2);
active_dt             = zeros(1,group_n);
FR_gr  = zeros(1,group_n);
CI95 = zeros(1,group_n);
for i = 1:group_n
    active_dt(1,i) = Analysis.groups{i}.active_per_s;
    tmp_neurons = Analysis.groups{i}.neurons;
    FR_gr(i) = mean(sum(sum(firings_store(:,tmp_neurons,:),1),3)/n_trial);
    CI95(i) = prctile(Perm.groups_perm{i},0.95);
end
active_dt_norm_FR_all  = active_dt/sum(FR);
active_dt_norm_FR_gr = active_dt./FR_gr;

CI95_norm_FR_all  = CI95/sum(FR);
CI95_norm_FR_gr = CI95./FR_gr;

if normalization == 0
    active_SO_dt     = active_dt;
    CI95_SO = (active_SO_dt >= CI95);
elseif normalization == 1
    active_SO_dt     = active_dt_norm_FR_all;
    CI95_SO = (active_SO_dt >= CI95_norm_FR_all);
elseif normalization == 2
    active_SO_dt     = active_dt_norm_FR_gr;
    CI95_SO = (active_SO_dt >= CI95_norm_FR_gr);
end
sum(CI95_SO)/group_n 

%%%%%%%%% activity vs acticity

mxx = max([active_W_dt active_SO_dt]);
bin_ind = 0:(mxx/50):mxx;

% activation W vs SO scatter
figure('Renderer', 'painters', 'Position', [10 10 260 500])
subplot(2,1,1)
scatter(active_W_dt, active_SO_dt,10,'.k')
hold on
plot([0 max(bin_ind)],[0 max(bin_ind)],'--k')
text(max(bin_ind)*0.75, max(bin_ind)*0.7, 'y=x');
xlabel('group activation (/s): Wake')
ylabel('group activation (/s): SO')
xlim([0 max(bin_ind)]); ylim([0 max(bin_ind)])

%%%%%%%%%%%%% histogram
subplot(2,1,2)
histogram(active_W_dt,bin_ind,'FaceColor',[0.6350 0.0780 0.1840])
hold on
histogram(active_SO_dt,bin_ind,'FaceColor',[0 0.4470 0.7410])
legend({'Wake','SO'})
xlabel('group activation (/s)')
ylabel('count');% set(gca, 'YScale', 'log')
xlim([0 max(bin_ind)])

%% group analysis; instantaneous activity W and SO

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
addpath('F:\20190624_Wisconsin\violin_plot')

% Analysis = load('20210118T020811_group_active_find_J0005.mat') % W
Analysis = load('20210118T021757_group_active_find_J00005.mat') % SO

firings_store    = Analysis.Data.trend.firings_store;
N                = Analysis.Data.param.N;
h                = Analysis.Data.param.h;
n_ind_s = round(1000/h);
n_trial = size(firings_store,1);

firing_long = zeros(N, n_ind_s*n_trial);
for s_i = 1:n_trial
    firing_long(:,((s_i-1)*n_ind_s+1):(s_i*n_ind_s)) = squeeze(firings_store(s_i,:,:));
end


kernel = gausswin(((2*12)/sqrt(2))/h);

len_act_inst = length(Analysis.groups{1}.active_inst);
t_ind = ((0:(len_act_inst-1))/1000)*h;
group_n   = size(Analysis.groups,2);
activity_inst_dt          = zeros(group_n,len_act_inst);

for i = 1:group_n
    
    disp([num2str(i)])
    activity_inst_dt(i,:)          = Analysis.groups{i}.active_inst;
end
activity_inst_dt_smooth         = filter(kernel,1,activity_inst_dt);
% 
% activity_inst_dt         = filter(kernel,1,activity_inst_dt);

activity_inst_all_dt         = sum(activity_inst_dt,1);

activity_inst_dt_bin = zeros(group_n,n_trial*1000);
for i = 1:(n_trial*1000 - 1)
    activity_inst_dt_bin(:,i) = sum(activity_inst_dt(:,((i-1)*(1/h)+1):((i)*(1/h))),2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% raster and inst active


plot_x_len = [0 1]; % SO

% plot_x_len = [7 8]; % W
plot_t_ind = min((plot_x_len*1000)/h + 1):max((plot_x_len*1000)/h + 1);

i_neur_1 = 1763; % 
i_neur_2 = 1680; % 
i_neur_3 = 2500; % 

figure
subplot(5,1,1)
[fire_I,fire_J] = find(firing_long);
fire_J = (fire_J*h)/1000;
scatter(fire_J,fire_I, 1,'.k');
% xlim([0 n_trial])
xlim(plot_x_len)
ylim([1 N])

subplot(5,1,2)
plot(t_ind,activity_inst_all_dt,'k')
% xlim([0 n_trial])
xlim(plot_x_len)
ylabel('iPR')

subplot(5,1,3)
% figure
i_neur = [i_neur_1];
% i_neur = [1300];
plot(t_ind,activity_inst_dt(i_neur,:),'k')
hold on
% hold on
% plot([0 n_trial],[0.5 0.5],'--k')
xlim(plot_x_len)
% ylim([0 1])
% ylabel({['group #' num2str(i_neur)], 'activation'},'r') %'Color',[0 0.4470 0.7410])
ylabel({'example group 1', 'activation'},'Color','r')
% legend('full group')
ylim([0 1])

subplot(5,1,4)
% figure
i_neur = [i_neur_2];
plot(t_ind,activity_inst_dt(i_neur,:),'k')
hold on
% hold on
% plot([0 n_trial],[0.5 0.5],'--k')
xlim(plot_x_len)
% ylim([0 1])
% ylabel({['group #' num2str(i_neur)], 'activation'},'b') %'Color',[0.6350 0.0780 0.1840])
ylabel({'example group 2', 'activation'},'Color','b')
% legend({'anchors','mid-exc','inh',['group #' num2str(i_neur)]})
ylim([0 1])

subplot(5,1,5)
% figure
i_neur = [i_neur_3];
plot(t_ind,activity_inst_dt(i_neur,:),'k')
hold on
% hold on
% plot([0 n_trial],[0.5 0.5],'--k')
xlim(plot_x_len)
% ylim([0 1])
% ylabel({['group #' num2str(i_neur)], 'activation'},'g') 'Color',[0.4660, 0.6740, 0.1880])
ylabel({'example group 3', 'activation'},'Color','g')
% legend({'anchors','mid-exc','inh',['group #' num2str(i_neur)]})
ylim([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% active vs active

% activity_inst_gr1_dt = activity_inst_dt_smooth(i_neur_1,:);
% activity_inst_gr2_dt = activity_inst_dt_smooth(i_neur_2,:);
% 
% figure 
% scatter(activity_inst_gr1_dt, activity_inst_gr2_dt,0.01,'.k')
% xlabel('group #1763 activity','Color',[0 0.4470 0.7410])
% ylabel('group #1680 activity','Color',[0.6350 0.0780 0.1840])
% xlim([0 1.4]); ylim([0 1.4]);
% 

%%%%%%%%%%%%%%%%%%%% PC
FC = NaN(group_n,group_n);
for i = 1:group_n
    disp([num2str(i)])
    for ii = (i+1):group_n
        disp([num2str(i) ' ' num2str(ii)])
        tmp_r  = corrcoef(activity_inst_dt_bin(i,:), activity_inst_dt_bin(ii,:));
%         tmp_r  = corrcoef(activity_inst_dt(i,(end/4):end), activity_inst_dt(ii,(end/4):end));
        FC(i,ii) = tmp_r(1,2);
    end
end
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_FC_W' ],'FC','-v7.3');
save([datestr(now, 'yyyymmddTHHMMSS') '_group_active_FC_SO' ],'FC','-v7.3');

figure
imagesc(FC)
caxis([0 1]) 
a = colorbar;
colormap(jet)
a.Label.String = 'FC';
xlabel('groups'); ylabel('groups')
%% time from group seq, type, PC, hi- low-FR analysis

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run

% wake
% Analysis = load('20201101T014517_group_active_find_J01.mat') % W
Analysis = load('20201224T182423_group_active_find_J0005.mat') % W
% Analysis = load('20201224T221541_group_active_find_J00005.mat') % SO

firings_store    = Analysis.Data.trend.firings_store;
N                = Analysis.Data.param.N;
h                = Analysis.Data.param.h;

n_ind_s = round(1000/h);
n_trial = size(firings_store,1);

firing_long = zeros(N, n_ind_s*n_trial);
for s_i = 1:n_trial
    firing_long(:,((s_i-1)*n_ind_s+1):(s_i*n_ind_s)) = squeeze(firings_store(s_i,:,:));
end

iPR = sum(firing_long,1);
kernel = gausswin(((2*12)/sqrt(2))/h);
iPR = filter(kernel,1,iPR);
iPR = iPR/N;

PC = zeros(1,N);
for i = 1:N
    r_tmp = corrcoef(iPR,firing_long(i,:));
    PC(i) = r_tmp(1,2);
end
PC_W = PC;

FR_W  = sum(firing_long,2)/n_trial;

% SO
% Analysis   = load('20201101T020003_group_active_find_J00004.mat') % SO
% % Analysis = load('20201224T182423_group_active_find_J0005.mat') % W
Analysis = load('20201224T221541_group_active_find_J00005.mat') % SO

firings_store    = Analysis.Data.trend.firings_store;
N                = Analysis.Data.param.N;
h                = Analysis.Data.param.h;

n_ind_s = round(1000/h);
n_trial = size(firings_store,1);

firing_long = zeros(N, n_ind_s*n_trial);
for s_i = 1:n_trial
    firing_long(:,((s_i-1)*n_ind_s+1):(s_i*n_ind_s)) = squeeze(firings_store(s_i,:,:));
end

iPR = sum(firing_long,1);
kernel = gausswin(((2*12)/sqrt(2))/h);
iPR = filter(kernel,1,iPR);
iPR = iPR/N;

PC = zeros(1,N);
for i = 1:N
    r_tmp = corrcoef(iPR,firing_long(i,:));
    PC(i) = r_tmp(1,2);
end
PC_SO = PC;

FR_SO  = sum(firing_long,2)/n_trial;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_N = Analysis.Data.param.A;
B_N = Analysis.Data.param.B;
C_N = Analysis.Data.param.C;

D_N = Analysis.Data.param.D;
group_n   = size(Analysis.groups,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fire_J_all = [];
fire_I_all = [];
for i = 1:group_n
    disp([num2str(i) ' '])
    firings = Analysis.groups{i}.firings;
    [fire_I,fire_J] = find(firings);
    fire_J = fire_J*h;
    
    fire_J_all = [fire_J_all; fire_J];
    fire_I_all = [fire_I_all; fire_I];
end

% I is neuron
% J is latency
% max(fire_J_all)dddddd

Avrg_fire = NaN(N,1);
for i = 1:N
    tmp_lat = fire_J_all(find(fire_I_all == i));
    tmp_lat(tmp_lat < 3) = [];
    Avrg_fire(i) = mean(tmp_lat);
end

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save('avrg_latency_group', 'Avrg_fire')
% figure
% histogram(Avrg_fire)

used_N = unique(fire_I_all);
fire_J_all_avrg = zeros(length(used_N),1);
fire_I_all_avrg = zeros(length(used_N),1);
for i = 1:length(used_N)
    fire_I_all_avrg(i) = used_N(i);
    fire_J_all_avrg(i) = mean(fire_J_all(find(fire_I_all == used_N(i)))) ;
end

%%% non-avrg
inh_ind = find(fire_I_all > N*0.8);
exc_ind = find(~(fire_I_all > N*0.8));

fire_J_all_exc = fire_J_all(exc_ind);
fire_J_all_inh = fire_J_all(inh_ind);

fire_I_all_exc = fire_I_all(exc_ind);
fire_I_all_inh = fire_I_all(inh_ind);


%%% average
inh_ind = find(fire_I_all_avrg > N*0.8);
exc_ind = find(~(fire_I_all_avrg > N*0.8));

fire_J_all_exc = fire_J_all_avrg(exc_ind);
fire_J_all_inh = fire_J_all_avrg(inh_ind);

fire_I_all_exc = fire_I_all_avrg(exc_ind);
fire_I_all_inh = fire_I_all_avrg(inh_ind);
%%%%%%%%%%%%%%%%%%%%%%%%%% within-group time, exc + inh

figure
subplot(1,6,1)
scatter(fire_J_all_exc, FR_W(fire_I_all_exc),20,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840])
hold on
scatter(fire_J_all_inh, FR_W(fire_I_all_inh),20,'.','MarkerEdgeColor',[0 0.4470 0.7410])
xlabel({'time fired in','group sequence(ms)'})
ylabel('Wake FR')
% legend({'exc','inh'})
ylim([0 80])
xlim([0 40])

subplot(1,6,2)
scatter(fire_J_all_exc, FR_SO(fire_I_all_exc),20,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840])
hold on
scatter(fire_J_all_inh, FR_SO(fire_I_all_inh),20,'.','MarkerEdgeColor',[0 0.4470 0.7410])
xlabel({'time fired in','group sequence(ms)'})
ylabel('SO FR')
ylim([0 80])
xlim([0 40])

subplot(1,6,3)
scatter(fire_J_all_exc, PC_W(fire_I_all_exc),20,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840])
hold on
scatter(fire_J_all_inh, PC_W(fire_I_all_inh),20,'.','MarkerEdgeColor',[0 0.4470 0.7410])
xlabel({'time fired in','group sequence(ms)'})
ylabel('Wake PC')
ylim([-0.02 0.14])
xlim([0 40])

subplot(1,6,4)
scatter(fire_J_all_exc, PC_SO(fire_I_all_exc),20,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840])
hold on
scatter(fire_J_all_inh, PC_SO(fire_I_all_inh),20,'.','MarkerEdgeColor',[0 0.4470 0.7410])
xlabel({'time fired in','group sequence(ms)'})
ylabel('SO PC')
ylim([-0.02 0.14])
xlim([0 40])

subplot(1,6,5)
scatter(fire_J_all_exc, D_N(fire_I_all_exc),20,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840])
xlabel({'time fired in','group sequence(ms)'})
ylabel('Izhikevich D-parameter')
% ylim([-0.02 0.14])
xlim([0 40])

subplot(1,6,6)
scatter(fire_J_all_inh, B_N(fire_I_all_inh),20,'.','MarkerEdgeColor',[0 0.4470 0.7410])
xlabel({'time fired in','group sequence(ms)'})
ylabel('Izhikevich B-parameter')
% ylim([-0.02 0.14])
xlim([0 40])

% figure
% histogram(fire_J_all_exc,0:0.05:25)
% set(gca, 'YScale', 'log')
% xlim([0 10])

%%%%%%%%%%%%%%%%%%%%% plot FR vs PC

figure
subplot(2,2,1)
scatter(FR_W, FR_SO,10,'.k')
xlabel('firing rate (/s) Wake')
ylabel('firing rate (/s) SO')
lsline

subplot(2,2,2)
scatter(PC_W, PC_SO,10,'.k')
xlabel('PC Wake')
ylabel('PC SO')
lsline

subplot(2,2,3)
scatter(FR_W, PC_W,10,'.k')
xlabel('firing rate (/s) Wake')
ylabel('PC Wake')
lsline

subplot(2,2,4)
scatter(FR_SO, PC_SO,10,'.k')
xlabel('firing rate (/s) SO')
ylabel('PC SO')
lsline

%%%%%%%%%%%%%%%%%%%%%%%%%%% neuron-type vs wake FR PC

figure
subplot(2,2,1)
scatter(D_N(1:800),FR_W(1:800), 20,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840])
ylabel('firing rate (/s) Wake')
xlabel('D-param')
lsline
ylim([0 80])

subplot(2,2,2)
scatter(B_N(801:1000),FR_W(801:1000), 20,'.','MarkerEdgeColor',[0 0.4470 0.7410])
ylabel('firing rate (/s) Wake')
xlabel('B-param')
lsline
ylim([0 80])

subplot(2,2,3)
scatter(D_N(1:800),PC_W(1:800), 20,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840])
ylabel('PC Wake')
xlabel('D-param')
lsline
ylim([-0.02 0.14])

subplot(2,2,4)
scatter(B_N(801:1000),PC_W(801:1000), 20,'.','MarkerEdgeColor',[0 0.4470 0.7410])
ylabel('PC Wake')
xlabel('B-param')
lsline
ylim([-0.02 0.14])

%%%%%%%%%%%%%%%%%%%%%%%%%%% neuron-type vs SO FR PC

figure
subplot(2,2,1)
scatter(D_N(1:800),FR_SO(1:800), 20,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840])
ylabel('firing rate (/s) SO')
xlabel('D-param')
lsline
ylim([0 10])

subplot(2,2,2)
scatter(B_N(801:1000),FR_SO(801:1000), 20,'.','MarkerEdgeColor',[0 0.4470 0.7410])
ylabel('firing rate (/s) SO')
xlabel('B-param')
lsline
ylim([0 10])

subplot(2,2,3)
scatter(D_N(1:800),PC_SO(1:800), 20,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840])
ylabel('PC SO')
xlabel('D-param')
lsline
ylim([0 0.11])

subplot(2,2,4)
scatter(B_N(801:1000),PC_SO(801:1000), 20,'.','MarkerEdgeColor',[0 0.4470 0.7410])
ylabel('PC SO')
xlabel('B-param')
lsline
ylim([0 0.11])



%% group analysis : neuron same connectivity

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
% load('20201103T215135_groups.mat') %
load('20201224T104121_groups_v2.mat') %

n_gr = length(groups);
gr_neurons = cell(1,n_gr);
mother_neuron = zeros(1,n_gr);
for i_gr = 1:n_gr
    disp(['gr ' num2str(i_gr)])
    i_firings = groups{i_gr}.firings;
    gr_neurons{i_gr} = find(sum(i_firings,2) > 0);
    
    mother_neuron(i_gr) = groups{i_gr}.mother_neuron;
end

conn_ne_same = zeros(n_gr, n_gr);
for i_gr = 1:n_gr
    disp(['gr ' num2str(i_gr)])
    for ii_gr = 1:n_gr
        
        mnn_neurons = min(length(gr_neurons{i_gr}),length(gr_neurons{ii_gr}));
        conn_ne_same(i_gr,ii_gr) = length(intersect(gr_neurons{i_gr},gr_neurons{ii_gr}))/mnn_neurons;
    end
end

figure
imagesc(conn_ne_same)
xlabel('groups')
ylabel('groups')
h = colorbar;
ylabel(h, 'portion of shared neurons')

unique_mother = unique(mother_neuron);
n_mother = length(unique_mother);

figure
imagesc(mother_neuron)
colormap(lines)

figure
plot(mother_neuron)
conn_ne_same_mother = zeros(n_mother, n_mother);
for i_gr = 1:n_mother
    disp(['mother ' num2str(i_gr)])
    for ii_gr = 1:n_mother
        
        i_uni = unique_mother(i_gr);
        ii_uni = unique_mother(ii_gr);
        i_mother_ind = find(mother_neuron == i_uni);
        ii_mother_ind = find(mother_neuron == ii_uni);
        
        conn_ne_same_mother(i_gr,ii_gr) = mean(mean(conn_ne_same(i_mother_ind,ii_mother_ind)));
    end
end

figure
imagesc(conn_ne_same_mother)

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_group_conn_ne_same_v2' ],'conn_ne_same','mother_neuron','-v7.3');
%

%% group analysis: mix, plot activity vs activity with color of mother neuron

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
addpath('F:\20190624_Wisconsin\violin_plot')
load('20201225T183751_group_conn_ne_same_v2.mat') % W

% wake
Analysis = load('20201224T182423_group_active_find_J0005.mat') % W

group_n   = size(Analysis.groups,2);
size_dt   = zeros(1,group_n);
path_dt   = zeros(1,group_n);
active_dt = zeros(1,group_n);
for i = 1:group_n
    size_dt(1,i)   = Analysis.groups{i}.size;
    path_dt(1,i)   = Analysis.groups{i}.longest_path;
    active_dt(1,i) = Analysis.groups{i}.active_per_s;
end

size_dt = size_dt;
path_dt = path_dt;
active_W_dt = active_dt;


Analysis = load('20201225T030332_group_active_find_J0005_FlipRast.mat') % W FlipRast

group_n   = size(Analysis.groups,2);
size_dt   = zeros(1,group_n);
path_dt   = zeros(1,group_n);
active_dt = zeros(1,group_n);
for i = 1:group_n
    size_dt(1,i)   = Analysis.groups{i}.size;
    path_dt(1,i)   = Analysis.groups{i}.longest_path;
    active_dt(1,i) = Analysis.groups{i}.active_per_s;
end

active_W_CT1_dt = active_dt;

Analysis = load('20201225T103235_group_active_find_J0005_FlipRast_ShuffN.mat') % W FlipRast ShuffN

group_n   = size(Analysis.groups,2);
size_dt   = zeros(1,group_n);
path_dt   = zeros(1,group_n);
active_dt = zeros(1,group_n);
for i = 1:group_n
    size_dt(1,i)   = Analysis.groups{i}.size;
    path_dt(1,i)   = Analysis.groups{i}.longest_path;
    active_dt(1,i) = Analysis.groups{i}.active_per_s;
end

active_W_CT2_dt = active_dt;

FR_W_mean = mean(Analysis.Data.trend.FR_trend);

% SO
Analysis   = load('20201224T221541_group_active_find_J00005.mat') % SO

group_n   = size(Analysis.groups,2);
size_dt   = zeros(1,group_n);
path_dt   = zeros(1,group_n);
active_dt = zeros(1,group_n);
for i = 1:group_n
    size_dt(1,i)   = Analysis.groups{i}.size;
    path_dt(1,i)   = Analysis.groups{i}.longest_path;
    active_dt(1,i) = Analysis.groups{i}.active_per_s;
end

active_SO_dt = active_dt;

Analysis   = load('20201225T065030_group_active_find_J00005_FlipRast.mat') % SO  FlipRast

group_n   = size(Analysis.groups,2);
size_dt   = zeros(1,group_n);
path_dt   = zeros(1,group_n);
active_dt = zeros(1,group_n);
for i = 1:group_n
    size_dt(1,i)   = Analysis.groups{i}.size;
    path_dt(1,i)   = Analysis.groups{i}.longest_path;
    active_dt(1,i) = Analysis.groups{i}.active_per_s;
end

active_SO_CT1_dt = active_dt;


Analysis   = load('20201225T141626_group_active_find_J00005_FlipRast_ShuffN.mat') % SO FlipRast ShuffN

group_n   = size(Analysis.groups,2);
size_dt   = zeros(1,group_n);
path_dt   = zeros(1,group_n);
active_dt = zeros(1,group_n);
for i = 1:group_n
    size_dt(1,i)   = Analysis.groups{i}.size;
    path_dt(1,i)   = Analysis.groups{i}.longest_path;
    active_dt(1,i) = Analysis.groups{i}.active_per_s;
end

active_SO_CT2_dt = active_dt;;

FR_SO_mean = mean(Analysis.Data.trend.FR_trend);
jit_ind = round(Analysis.Data.param.gr_active_jit/Analysis.Data.param.h)*2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot
active_W_dt_norm      = active_W_dt     /(FR_W_mean*jit_ind);
active_W_CT1_dt_norm  = active_W_CT1_dt /(FR_W_mean*jit_ind);
active_W_CT2_dt_norm  = active_W_CT2_dt /(FR_W_mean*jit_ind);
active_SO_dt_norm     = active_SO_dt    /(FR_SO_mean*jit_ind);
active_SO_CT1_dt_norm = active_SO_CT1_dt/(FR_SO_mean*jit_ind);
active_SO_CT2_dt_norm = active_SO_CT2_dt/(FR_SO_mean*jit_ind);
active_diff_dt_norm   = active_SO_dt_norm - active_W_dt_norm;

active_W_dt_norm_sz = [];
active_SO_dt_norm_sz = [];
active_diff_dt_norm_sz = [];
for sz = min(size_dt):max(size_dt)
    active_W_dt_norm_sz{sz-min(size_dt)+1} =active_W_dt_norm(find(size_dt == sz));
    active_SO_dt_norm_sz{sz-min(size_dt)+1} =active_SO_dt_norm(find(size_dt == sz));
    active_diff_dt_norm_sz{sz-min(size_dt)+1} =active_diff_dt_norm(find(size_dt == sz));
end

%%%%%%%%% activity vs acticity

unique_mother = unique(mother_neuron);
n_mother = length(unique_mother);

figure
subplot(1,3,1)
% scatter(active_W_dt_norm, active_SO_dt_norm,10,'.k')
% hold on
for ii = 1:n_mother
    mother_ii = find(mother_neuron == unique_mother(ii));
    scatter(active_W_dt_norm(mother_ii), active_SO_dt_norm(mother_ii),20,'.')
    hold on
end
plot([0 10],[0 10],'--k')
text(3.6, 3.3, 'y=x');
xlabel('group activation (/s): Wake')
ylabel('group activation (/s): SO')
xlim([0 8])
ylim([0 3.5])
title('observed')

subplot(1,3,2)
% scatter(active_W_CT1_dt_norm, active_SO_CT1_dt_norm,10,'.k')
% hold on
for ii = 1:n_mother
    mother_ii = find(mother_neuron == unique_mother(ii));
    scatter(active_W_CT1_dt_norm(mother_ii), active_SO_CT1_dt_norm(mother_ii),20,'.')
    hold on
end
plot([0 10],[0 10],'--k')
% text(3.6, 3.3, 'y=x');
xlabel('group activation (/s): Wake')
ylabel('group activation (/s): SO')
xlim([0 8])
ylim([0 3.5])
title('control: flip raster ')

subplot(1,3,3)
for ii = 1:n_mother
    mother_ii = find(mother_neuron == unique_mother(ii));
    scatter(active_W_CT2_dt_norm(mother_ii), active_SO_CT2_dt_norm(mother_ii),20,'.')
    hold on
end
plot([0 10],[0 10],'--k')
% text(3.6, 3.3, 'y=x');
xlabel('group activation (/s): Wake')
ylabel('group activation (/s): SO')
xlim([0 8])
ylim([0 3.5])
title('control: flip raster + shuffle neurons')

%% group analysis: group type mother, anchor
%%%%%%%%%% need update J0005, new group v2
clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
addpath('F:\20190624_Wisconsin\violin_plot')
load('20201104T162112_group_conn_ne_same.mat')

% wake
Analysis = load('20201103T215135_groups.mat') % % W

group_n   = size(Analysis.groups,2);
size_dt   = zeros(1,group_n);
path_dt   = zeros(1,group_n);
mother_dt   = zeros(1,group_n);
mother_used_dt   = zeros(1,group_n);
anchor_dt = zeros(3,group_n);
% non_anchor_non_mother_ex_dt = cell(1,group_n);
gr_neurons_dt = cell(1,group_n);
for i = 1:group_n
    size_dt(1,i)   = Analysis.groups{i}.size;
    path_dt(1,i)   = Analysis.groups{i}.longest_path;
    anchor_dt(:,i) = Analysis.groups{i}.anchor_neurons;
    mother_dt(1,i) = Analysis.groups{i}.mother_neuron;
    
    i_firings = Analysis.groups{i}.firings;
    gr_neurons_dt{i} = find(sum(i_firings,2) > 0);
    
    mother_used_dt(1,i) = length(find(gr_neurons_dt{i} == mother_dt(1,i))) > 0;
    
    %     tmp_neur = setdiff(gr_neurons_dt{i},anchor_dt(:,i));
    %     tmp_neur = setdiff(tmp_neur,mother_dt(1,i));
    %     non_anchor_non_mother_ex_dt{1,i} = tmp_neur(find(tmp_neur < 801)); %%%%%%%% removes repeated neru too
end

sum(mother_used_dt)/group_n

figure
imagesc(mother_used_dt)
colormap(gray)

size_dt = size_dt;
path_dt = path_dt;
active_W_dt = active_dt;
%%%%%%%%%%%%%%%%%%%%%%%%

C_tmp = Analysis.param.C;
D_tmp = Analysis.param.D;

unique_mother = unique(mother_neuron);
n_mother = length(unique_mother);

unique_anchor = unique(anchor_dt);

figure
histogram(D_tmp(1:800),2:1:8,'FaceColor',[0.7 0.7 0.7])
hold on
histogram(D_tmp(unique_mother),2:1:8,'FaceColor',[0.6350 0.0780 0.1840])
xlabel('Izhikevich parameter d')
ylabel('neuron count')
legend({'excitatory neurons','mother neurons'},'Location','Best')
xlim([2 8])

figure
histogram(D_tmp(1:800),2:1:8,'FaceColor',[0.7 0.7 0.7])
hold on
histogram(D_tmp(unique_anchor),2:1:8,'FaceColor',[0.6350 0.0780 0.1840])
xlabel('Izhikevich parameter d')
ylabel('neuron count')
legend({'excitatory neurons','anchor neurons'},'Location','Best')
xlim([2 8])

%% stepwise decrease states 100s, 0.1 to 0.0005

clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
load('20201028T133420_Izhi_run.mat') %beginning of wake
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
addpath('E:\poissonTutorial')

J_ATP_step = [0.1 0.05 0.01 0.005 0.001 0.0005];

for i = 1:length(J_ATP_step)
    
    disp(['J ' num2str(J_ATP_step(i))])
    
    initial.T_sec   = 100; %sec
    param.stdp1     = 0;
    param.display   = 0;
    
    param.J_ATP   = J_ATP_step(i);
    
    [initial, trend] = Izhi_ATP_run_v6(initial,param);
    
    cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
    save([datestr(now, 'yyyymmddTHHMMSS') '_Izhi_run100s_J_step' num2str(i) ],'-v7.3');
    
end

%% stepwise decrease states 100s, 0.005 to 0.0005

clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
load('20201107T023748_Izhi_run100s_J_step4.mat') % after 100s J = 0.005
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
addpath('E:\poissonTutorial')

J_ATP_step = [0.005 0.0043 0.0035 0.0028 0.002 0.0013 0.0005];

for i = 1:length(J_ATP_step)
    
    disp(['J ' num2str(J_ATP_step(i))])
    
    initial.T_sec   = 100; %sec
    param.stdp1     = 0;
    param.display   = 0;
    
    param.J_ATP   = J_ATP_step(i);
    
    [initial, trend] = Izhi_ATP_run_v6(initial,param);
    
    cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
    save([datestr(now, 'yyyymmddTHHMMSS') '_Izhi_run100s_J0005to00005_step' num2str(i) ],'-v7.3');
    
end


%% stepwise decrease states 100s plot, 0.005 to 0.0005



clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
load('20201028T133420_Izhi_run.mat') %beginning of wake
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
addpath('E:\poissonTutorial')

% J_ATP_step = [0.1 0.05 0.01 0.005 0.001 0.0005];
J_ATP_step = [0.005 0.0043 0.0035 0.0028 0.002 0.0013 0.0005];
n_step = length(J_ATP_step);

C_long = zeros(n_step, (1000/0.1)*100);

figure
for i = 1:n_step
    fileList = dir(['*_Izhi_run100s_J0005to00005_step' num2str(i) '.mat']);
    %     fileList = dir(['*_Izhi_run100s_J_step' num2str(i) '.mat']);
    load(fileList.name);
    
    firings_store    = trend.firings_store;
    N                = param.N;
    h                = param.h;
    
    n_ind_s = round(1000/h);
    n_trial = size(firings_store,1);
    
    firing_long = zeros(N, n_ind_s*n_trial);
    for s_i = 1:n_trial
        firing_long(:,((s_i-1)*n_ind_s+1):(s_i*n_ind_s)) = squeeze(firings_store(s_i,:,:));
    end
    
    tic
    subplot(n_step,1,i)
    firing_long = zeros(N, n_ind_s*n_trial);
    for s_i = 1:n_trial
        firing_long(:,((s_i-1)*n_ind_s+1):(s_i*n_ind_s)) = squeeze(firings_store(s_i,:,:));
    end
    [fire_I,fire_J] = find(firing_long);
    fire_J = (fire_J*h)/1000;
    scatter(fire_J,fire_I, 1,'.k');
    xlim([0 n_trial])
    %     xlim([80 n_trial])
    ylim([1 N])
    toc
    
    C_ATP_store_1N    = trend.C_ATP_store_1N;
    
    for s_i = 1:n_trial
        C_long(i,((s_i-1)*n_ind_s+1):(s_i*n_ind_s)) = C_ATP_store_1N(s_i,:);
    end
    
    
end

mnn_C = min(min(C_long));
mxx_C = max(max(C_long));
figure
for i = 1:n_step
    C_long_i = C_long(i,:);
    
    t_ind = ((1:length(C_long))*h)/1000;
    
    subplot(n_step,1,i)
    plot(t_ind, C_long_i,'LineWidth',1)
    ylim([mnn_C mxx_C])
    
end

figure
plot([1 repelem(2:n_step,2) n_step+1], repelem(J_ATP_step,2),'k')
ylim([0 0.006])
% set(gca, 'YScale', 'log')
xticks(1:(n_step+1))
xticklabels(0:100:600)
xlabel('s'); ylabel('J_{ATP}')

%% Wake: stimuli J 0.005 (rand stim)

clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
%%%%%%%%%%%%%%%%%%%%
load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
% load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
addpath('E:\poissonTutorial')

% load('20201228T204604_letter_stim_neurons_randi.mat') % instance of stimuli, letter

initial.T_sec   = 100; %sec
param.stdp1     = 0;
param.display   = 1;

param.J_ATP   = 0.005; %
% param.J_ATP   = 0.0005; %

stim_neurons         = randi(1000,[30 1]);
letter               = zeros(1000,1);
letter(stim_neurons) = 1;
save([datestr(now, 'yyyymmddTHHMMSS') '_letter_stim_neurons_randi' ],'letter','stim_neurons');
% 20201228T204604_letter_stim_neurons_randi

initial.letter       = letter;
initial.stim_neurons = stim_neurons;

param.stimuli_onset_1s                        = zeros(1,1000/(param.h));
param.stimuli_onset_1s([250 750]/(param.h))   = 1;

[initial, trend] = Izhi_ATP_run_v6(initial,param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_Izhi_run_J0005_stimuli_randi_v4' ],'-v7.3');
% 20201228T235903_Izhi_run_J0005_stimuli_randi_v3
% 20201231T192605_Izhi_run_J0005_stimuli_randi_v4

%% SO: stimuli J 0.0005 (rand stim)

clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
%%%%%%%%%%%%%%%%%%%%
% load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
addpath('E:\poissonTutorial')

load('20201228T204604_letter_stim_neurons_randi.mat') % instance of stimuli, letter

initial.T_sec   = 100; %sec
param.stdp1     = 0;
param.display   = 1;

% param.J_ATP   = 0.005; %
param.J_ATP   = 0.0005; %

initial.letter       = letter;
initial.stim_neurons = stim_neurons;

param.stimuli_onset_1s                        = zeros(1,1000/(param.h));
param.stimuli_onset_1s([250 750]/(param.h))   = 1;

[initial, trend] = Izhi_ATP_run_v6(initial,param);

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_Izhi_run_J00005_stimuli_randi_v4' ],'-v7.3');
% 20201229T031229_Izhi_run_J00005_stimuli_randi_v3
% 20201231T225045_Izhi_run_J00005_stimuli_randi_v4


%% PETH, stimuli Wake J 0.005 vs. SO J 0.0005

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
    excl_len      = 10/h; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAM
    
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

%% PETH, spontaneous , SO J 0.0005, UP detection

clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')

% load('20201108T224354_Izhi_run100s_J0005to00005_step1.mat') % after 100s J 0.005
Analysis   = load('20201109T043258_Izhi_run100s_J0005to00005_step7.mat') % after 100s J 0.0005

firings_store    = Analysis.trend.firings_store;
N                = Analysis.param.N;
h                = Analysis.param.h;

lag_segment = 0/h;
pre_segment = 0/h;
segment = 200/h;

convFilter =gausswin(4);

n_ind_s = round(1000/h);
n_trial = size(firings_store,1);

firing_long = zeros(N, n_ind_s*n_trial);
for s_i = 1:n_trial
    firing_long(:,((s_i-1)*n_ind_s+1):(s_i*n_ind_s)) = squeeze(firings_store(s_i,:,:));
end

iPR = sum(firing_long,1);
kernel = gausswin(((2*12)/sqrt(2))/h);
iPR = filter(kernel,1,iPR);
iPR = iPR/N;

PC = zeros(1,N);
for i = 1:N
    r_tmp = corrcoef(iPR,firing_long(i,:));
    PC(i) = r_tmp(1,2);
end
PC_SO = PC;

FR_SO  = sum(firing_long,2)/n_trial;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_N = Analysis.param.A;
B_N = Analysis.param.B;
C_N = Analysis.param.C;

D_N = Analysis.param.D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Spk_train = firing_long';

spk_sum = sum(Spk_train,2);
UP_60ms_min = floor(size(Spk_train,2)*(0.4));
%          UP_60ms_min = 15;
data_len = size(Spk_train,1);

spk_sum_60ms = conv(ones(60/h,1),spk_sum);
spk_sum_60ms = spk_sum_60ms((length(spk_sum_60ms) - data_len + 1):end);
spk_sum_30ms = conv(spk_sum,ones(30/h,1));
spk_sum_30ms = spk_sum_30ms(1:(end - (length(spk_sum_30ms) - data_len)));
tmp_UP         = (spk_sum_30ms <= 1).*(spk_sum_60ms >= UP_60ms_min);
tmp_spk        = find(tmp_UP)';
tmp_x          = [diff(tmp_spk)~=1,true];
ind_UP         = tmp_spk(tmp_x);
ind_UP(find(diff(ind_UP) < 100)+1) = [];

ind_UP(ind_UP <= pre_segment) = [];
ind_UP(ind_UP >= (data_len-segment))     = [];
spk_sum_UP         = zeros(size(spk_sum));
spk_sum_UP(ind_UP) = 1;

bin = 1/h;

Spk_UP = zeros(length(ind_UP),segment+pre_segment-lag_segment+bin,size(Spk_train,2));
for i = 1:length(ind_UP)
    Spk_UP(i,:,:) = Spk_train((ind_UP(i)-pre_segment+lag_segment):(ind_UP(i)+segment-1+bin),:);
end


Spk_UP_bin = zeros(size(Spk_UP,1),((size(Spk_UP,2)-bin)/bin),size(Spk_UP,3));
for i = 1:((size(Spk_UP,2)-bin)/bin)
    Spk_UP_bin(:,i,:) = sum(Spk_UP(:,((i-1)*bin + 1):(i*bin),:),2);
end

% figure
% tmp = squeeze(sum(Spk_UP_bin,1));
% imagesc(tmp')


Spk_UP_perm = permute(Spk_UP_bin,[2 3 1]);


% PETH
Spk_train_conv_UP = [];
for i = 1:size(Spk_UP_perm,2)
    Spk_train_conv_UP(:,i) = conv(sum(Spk_UP_perm(:,i,:),3),convFilter);
end


tmp_mxx = max(Spk_train_conv_UP,[],1);
tmp_mnn = min(Spk_train_conv_UP,[],1);
%         tmp_mxx = max(Spk_train_conv_UP(:));
%         tmp_mnn = min(Spk_train_conv_UP(:));
Spk_train_conv_UP_norm = (Spk_train_conv_UP - tmp_mnn)./(tmp_mxx - tmp_mnn);


cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save(['20201220_Spk_UP_J00005'],'Spk_UP_bin')
save(['20201220_FR_SO_J00005'],'FR_SO')
save(['20201220_PC_SO_J00005'],'PC_SO')

%%%%%%%%%%%%%%%%%%%%%%
%
% figure
% imagesc(Spk_train_conv_UP_norm')
%
%
% %%%%%%%%%%%%%%
% ind = 35000:45000;
%
% figure
% subplot(5,1,1)
% imagesc(Spk_train(ind,:)')
% ylabel('neurons')
%
% subplot(5,1,2)
% plot(spk_sum(ind))
% xlim([0 length(ind)])
% ylabel('iPR')
%
% subplot(5,1,3)
% plot(spk_sum_30ms(ind))
% hold on
% plot((spk_sum_30ms(ind) <= 1)*max(spk_sum_30ms(ind)),'--k','LineWidth',1)
% xlim([0 length(ind)])
% ylabel('conv 30ms before')
%
% subplot(5,1,4)
% plot(spk_sum_60ms(ind))
% hold on
% plot((spk_sum_60ms(ind) >= UP_60ms_min)*max(spk_sum_60ms(ind)),'--k','LineWidth',1)
% xlim([0 length(ind)])
% ylabel('conv 60ms after')
%
% subplot(5,1,5)
% plot(spk_sum_UP(ind),'k')
% xlim([0 length(ind)])
% ylabel('UP onset')
% xlabel('ms')
% % subplot(5,1,5)
% % plot(tmp_UP(ind),'k')
% % xlim([0 length(ind)])

%% PETH: compare UP with stiuli W SO

clear all

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run

% cenmass_SO                = load('20201231_Izhi_run_J00005_SO_stimuli_cenmass.mat').cenmass;
% firings_store_prestim_SO  = load('20201231_Izhi_run_J00005_SO_stimuli_firings_store_prestim.mat').firings_store_prestim;
% firings_store_poststim_SO = load('20201231_Izhi_run_J00005_SO_stimuli_firings_store_poststim.mat').firings_store_poststim;
% sig_neuron_SO             = load('20201231_Izhi_run_J00005_SO_stimuli_sig_neuron.mat').sig_neuron;
% 
% cenmass_W                 = load('20201231_Izhi_run_J0005_W_stimuli_cenmass.mat').cenmass;
% firings_store_poststim_W  = load('20201231_Izhi_run_J0005_W_stimuli_firings_store_poststim.mat').firings_store_poststim;
% sig_neuron_W              = load('20201231_Izhi_run_J0005_W_stimuli_sig_neuron.mat').sig_neuron;
% 
% data_trans                = load('20201220_Spk_UP_J00005_trans.mat').data_trans;

cenmass_SO                = load('20201221_Izhi_run_J00005_SO_stimuli_cenmass.mat').cenmass;
firings_store_prestim_SO  = load('20201221_Izhi_run_J00005_SO_stimuli_firings_store_prestim.mat').firings_store_prestim;
firings_store_poststim_SO = load('20201221_Izhi_run_J00005_SO_stimuli_firings_store_poststim.mat').firings_store_poststim;
sig_neuron_SO             = load('20201221_Izhi_run_J00005_SO_stimuli_sig_neuron.mat').sig_neuron;

cenmass_W                 = load('20201221_Izhi_run_J0005_W_stimuli_cenmass.mat').cenmass;
firings_store_poststim_W  = load('20201221_Izhi_run_J0005_W_stimuli_firings_store_poststim.mat').firings_store_poststim;
sig_neuron_W              = load('20201221_Izhi_run_J0005_W_stimuli_sig_neuron.mat').sig_neuron;

data_trans                = load('20201220_Spk_UP_J00005_trans.mat').data_trans;


convFilter =gausswin(4);

Spk_UP_perm = permute(data_trans,[2 3 1]);

% PETH
Spk_train_conv_UP = [];
for i = 1:size(Spk_UP_perm,2)
    Spk_train_conv_UP(:,i) = conv(sum(Spk_UP_perm(:,i,:),3),convFilter);
end

excl_len      = 10; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAM
N = size(Spk_train_conv_UP,2);
conv_sig = (Spk_train_conv_UP > prctile(firings_store_prestim_SO',95))'; %%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAM
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


%%%%%% need to redo later, match methods
bin = 1/0.1;
Spk_poststim_W_bin = zeros(N,100);
Spk_poststim_SO_bin = zeros(N,100);
for i = 1:99
    Spk_poststim_W_bin(:,i) = sum(firings_store_poststim_W(:,((i-1)*bin+1):(i*bin)),2);
    Spk_poststim_SO_bin(:,i) = sum(firings_store_poststim_SO(:,((i-1)*bin+1):(i*bin)),2);
end


tmp_mxx = max(Spk_train_conv_UP,[],1);
tmp_mnn = min(Spk_train_conv_UP,[],1);
Spk_train_conv_UP_norm = (Spk_train_conv_UP - tmp_mnn)./(tmp_mxx - tmp_mnn);

tmp_mxx = max(Spk_poststim_W_bin,[],1);
tmp_mnn = min(Spk_poststim_W_bin,[],1);
Spk_poststim_W_norm = (Spk_poststim_W_bin - tmp_mnn)./(tmp_mxx - tmp_mnn);

tmp_mxx = max(Spk_poststim_SO_bin,[],1);
tmp_mnn = min(Spk_poststim_SO_bin,[],1);
Spk_poststim_SO_norm = (Spk_poststim_SO_bin - tmp_mnn)./(tmp_mxx - tmp_mnn);


conv_sig_firing = conv_sig'.*Spk_train_conv_UP_norm((1:end-1),:);

segment = 100;
cenmass_UP = [];
for i = 1:size(conv_sig_firing,2)
    cenmass_UP(i) = sum(conv_sig_firing(1:segment,i)'.*(1:segment))/sum(conv_sig_firing(1:segment,i));
end
sig_neuron_UP = sum(conv_sig,2) ~= 0 ;
cenmass_UP_orig = cenmass_UP;
cenmass_W_orig = cenmass_W;
cenmass_SO_orig = cenmass_SO;
%
sig_neuron_all = ((sig_neuron_SO + sig_neuron_W + sig_neuron_UP) == 3);
% 
[~,Ind_sort] = sort(cenmass_W);
% [~,Ind_sort] = sort(cenmass_UP);
Spk_poststim_W_norm    = Spk_poststim_W_norm(Ind_sort,:);
Spk_poststim_SO_norm   = Spk_poststim_SO_norm(Ind_sort,:);
Spk_train_conv_UP_norm = Spk_train_conv_UP_norm(1:100,Ind_sort)';
cenmass_W      = cenmass_W(Ind_sort); %/bin;
cenmass_SO     = cenmass_SO(Ind_sort); %/bin;
cenmass_UP     = cenmass_UP(Ind_sort)';
sig_neuron_all = sig_neuron_all(Ind_sort);
sig_neuron_ind = find(sig_neuron_all);

I_EI = [ones(800,1); zeros(200,1)];
I_EI = I_EI(Ind_sort);
Ind_E = find(I_EI.*sig_neuron_all);
Ind_I = find((~I_EI).*sig_neuron_all);

Ind_all = 1:length(sig_neuron_ind);
Ind_E_pos = [];
Ind_I_pos = [];
for i = 1:length(sig_neuron_ind);
    if any(Ind_E == sig_neuron_ind(i))
        Ind_E_pos = [Ind_E_pos i];
    elseif any(Ind_I == sig_neuron_ind(i))
        Ind_I_pos = [Ind_I_pos i];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,3,1)
imagesc(Spk_poststim_W_norm(sig_neuron_ind,:))
hold on
scatter(cenmass_W(sig_neuron_ind), 1:length(sig_neuron_ind),30,'+w')
xlabel('time from stimuli (ms)'); ylabel('neurons')
title('Wake(J=0.005):stimuli-response')
caxis([0 1])

subplot(2,3,2)
imagesc(Spk_poststim_SO_norm(sig_neuron_ind,:))
hold on
scatter(cenmass_SO(sig_neuron_ind), 1:length(sig_neuron_ind),30,'+w')
xlabel('time from stimuli (ms)'); 
title('SO(J=0.0005):stimuli-response')
caxis([0 1])

subplot(2,3,3)
imagesc(Spk_train_conv_UP_norm(sig_neuron_ind,:))
hold on
scatter(cenmass_UP(sig_neuron_ind), 1:length(sig_neuron_ind),30,'+w')
xlabel('time from UP-onset (ms)'); 
title('SO(J=0.0005):UP-state')
caxis([0 1])

subplot(2,3,4)
scatter(cenmass_W(Ind_E), Ind_E_pos,100,'.r')
hold on
scatter(cenmass_W(Ind_I), Ind_I_pos,100,'.b')
set(gca, 'YDir','reverse');xlim([0 100])
caxis([0 1])

subplot(2,3,5)
scatter(cenmass_SO(Ind_E), Ind_E_pos,100,'.r')
hold on
scatter(cenmass_SO(Ind_I), Ind_I_pos,100,'.b')
set(gca, 'YDir','reverse');xlim([0 100])
caxis([0 1])

subplot(2,3,6)
scatter(cenmass_UP(Ind_E), Ind_E_pos,100,'.r')
hold on
scatter(cenmass_UP(Ind_I), Ind_I_pos,100,'.b')
set(gca, 'YDir','reverse');xlim([0 100])
caxis([0 1])

%%%%%%%%%%%%%%%% SD
figure
subplot(1,3,1)
scatter(cenmass_SO(Ind_E) - cenmass_W(Ind_E), Ind_E_pos,100,'.r')
hold on
scatter(cenmass_SO(Ind_I) - cenmass_W(Ind_I), Ind_I_pos,100,'.b')
SD_E = std(cenmass_SO(Ind_E) - cenmass_W(Ind_E))
SD_I = std(cenmass_SO(Ind_I) - cenmass_W(Ind_I))
title(['\color{red}SD_{exc}=' num2str(SD_E) ', \color{blue}SD_{inh}=' num2str(SD_I)])
xlim([-100 100]);  ylim([1 length(sig_neuron_ind)])
ylabel('neurons')
xlabel('Latency_{SO,stimuli} - Latency_{Wake,stimuli}(ms)')

subplot(1,3,2)
scatter(cenmass_UP(Ind_E) - cenmass_W(Ind_E), Ind_E_pos,100,'.r')
hold on
scatter(cenmass_UP(Ind_I) - cenmass_W(Ind_I), Ind_I_pos,100,'.b')
SD_E = std(cenmass_UP(Ind_E) - cenmass_W(Ind_E))
SD_I = std(cenmass_UP(Ind_I) - cenmass_W(Ind_I))
title(['\color{red}SD_{exc}=' num2str(SD_E) ', \color{blue}SD_{inh}=' num2str(SD_I)])
xlim([-100 100]); ylim([1 length(sig_neuron_ind)])
xlabel('Latency_{SO,UP} - Latency_{Wake,stimuli}(ms)')

subplot(1,3,3)
scatter(cenmass_SO(Ind_E) - cenmass_UP(Ind_E), Ind_E_pos,100,'.r')
hold on
scatter(cenmass_SO(Ind_I) - cenmass_UP(Ind_I), Ind_I_pos,100,'.b')
SD_E = std(cenmass_SO(Ind_E) - cenmass_UP(Ind_E))
SD_I = std(cenmass_SO(Ind_I) - cenmass_UP(Ind_I))
title(['\color{red}SD_{exc}=' num2str(SD_E) ', \color{blue}SD_{inh}=' num2str(SD_I)])
xlim([-100 100]); ylim([1 length(sig_neuron_ind)])
ylabel('neurons')
xlabel('Latency_{SO,stimuli} - Latency_{SO,UP}(ms)')

%%%%%%%%%% corr UPlatency with avrg group latency

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
load('avrg_latency_group.mat')

cen_orig = cenmass_SO_orig;
% cen_orig = cenmass_UP_orig';

not_nan = find(~isnan(Avrg_fire).*~isnan(cen_orig));

figure
scatter(Avrg_fire(not_nan), cen_orig(not_nan))




%% PETH: group activity: stimuli-W, stimuli-SO, UP-SO

clear all
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run

Groups = load('20201224T104121_groups_v2.mat') %

cenmass_SO                = load('20201231_Izhi_run_J00005_SO_stimuli_cenmass.mat').cenmass;
firings_store_prestim_SO  = load('20201231_Izhi_run_J00005_SO_stimuli_firings_store_prestim.mat').firings_store_prestim;
firings_store_poststim_SO = load('20201231_Izhi_run_J00005_SO_stimuli_firings_store_poststim.mat').firings_store_poststim;
sig_neuron_SO             = load('20201231_Izhi_run_J00005_SO_stimuli_sig_neuron.mat').sig_neuron;
% sig_neuron_SO             = load('20201221_Izhi_run_J00005_SO_stimuli_sig_neuron.mat').sig_neuron;

cenmass_W                 = load('20201231_Izhi_run_J0005_W_stimuli_cenmass.mat').cenmass;
firings_store_poststim_W  = load('20201231_Izhi_run_J0005_W_stimuli_firings_store_poststim.mat').firings_store_poststim;
sig_neuron_W              = load('20201231_Izhi_run_J0005_W_stimuli_sig_neuron.mat').sig_neuron;
% sig_neuron_W              = load('20201221_Izhi_run_J0005_W_stimuli_sig_neuron.mat').sig_neuron;

data_trans                = load('20201220_Spk_UP_J00005_trans.mat').data_trans;


convFilter =gausswin(4);

Spk_UP_perm = permute(data_trans,[2 3 1]);

% PETH
Spk_train_conv_UP = [];
for i = 1:size(Spk_UP_perm,2)
    Spk_train_conv_UP(:,i) = conv(sum(Spk_UP_perm(:,i,:),3),convFilter);
end

excl_len      = 10;
N = size(Spk_train_conv_UP,2);
conv_sig = (Spk_train_conv_UP > prctile(firings_store_prestim_SO',95))';
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


%%%%%% need to redo later, match methods
bin = 1/0.1;
Spk_poststim_W_bin = zeros(N,100);
Spk_poststim_SO_bin = zeros(N,100);
for i = 1:99
    Spk_poststim_W_bin(:,i) = sum(firings_store_poststim_W(:,((i-1)*bin+1):(i*bin)),2);
    Spk_poststim_SO_bin(:,i) = sum(firings_store_poststim_SO(:,((i-1)*bin+1):(i*bin)),2);
end


tmp_mxx = max(Spk_train_conv_UP,[],1);
tmp_mnn = min(Spk_train_conv_UP,[],1);
Spk_train_conv_UP_norm = (Spk_train_conv_UP - tmp_mnn)./(tmp_mxx - tmp_mnn);

tmp_mxx = max(Spk_poststim_W_bin,[],1);
tmp_mnn = min(Spk_poststim_W_bin,[],1);
Spk_poststim_W_norm = (Spk_poststim_W_bin - tmp_mnn)./(tmp_mxx - tmp_mnn);

tmp_mxx = max(Spk_poststim_SO_bin,[],1);
tmp_mnn = min(Spk_poststim_SO_bin,[],1);
Spk_poststim_SO_norm = (Spk_poststim_SO_bin - tmp_mnn)./(tmp_mxx - tmp_mnn);


conv_sig_firing = conv_sig'.*Spk_train_conv_UP_norm((1:end-1),:);

segment = 100;
cenmass_UP = [];
for i = 1:size(conv_sig_firing,2)
    cenmass_UP(i) = sum(conv_sig_firing(1:segment,i)'.*(1:segment))/sum(conv_sig_firing(1:segment,i));
end
sig_neuron_UP = sum(conv_sig,2) ~= 0 ;
%
sig_neuron_all = ((sig_neuron_SO + sig_neuron_W + sig_neuron_UP) == 3);



[~,Ind_sort] = sort(cenmass_W);
Spk_poststim_W_norm    = Spk_poststim_W_norm(Ind_sort,:);
Spk_poststim_SO_norm   = Spk_poststim_SO_norm(Ind_sort,:);
Spk_train_conv_UP_norm = Spk_train_conv_UP_norm(1:100,Ind_sort)';
cenmass_W      = cenmass_W(Ind_sort); %/bin;
cenmass_SO     = cenmass_SO(Ind_sort); %/bin;
cenmass_UP     = cenmass_UP(Ind_sort)';

sig_neuron_all_orig = sig_neuron_all;
sig_neuron_ind_orig = find(sig_neuron_all_orig);

sig_neuron_all = sig_neuron_all(Ind_sort);
sig_neuron_ind = find(sig_neuron_all);

I_EI = [ones(800,1); zeros(200,1)];
I_EI = I_EI(Ind_sort);
Ind_E = find(I_EI.*sig_neuron_all);
Ind_I = find((~I_EI).*sig_neuron_all);

Ind_all = 1:length(sig_neuron_ind);
Ind_E_pos = [];
Ind_I_pos = [];
for i = 1:length(sig_neuron_ind);
    if any(Ind_E == sig_neuron_ind(i))
        Ind_E_pos = [Ind_E_pos i];
    elseif any(Ind_I == sig_neuron_ind(i))
        Ind_I_pos = [Ind_I_pos i];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% 

groups              = Groups.groups;
param               = Groups.param;
param.gr_active_jit = 1;

peth = Spk_poststim_SO_norm;
% peth = Spk_poststim_W_norm;
% peth = Spk_train_conv_UP_norm;

[groups_active] = group_active_find_PETH_v0(groups,param,peth,sig_neuron_ind_orig)
% groups_active: neurons(gr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,3,1)
imagesc(Spk_poststim_W_norm(sig_neuron_ind,:))
hold on
scatter(cenmass_W(sig_neuron_ind), 1:length(sig_neuron_ind),30,'+w')
xlabel('time from stimuli (ms)'); ylabel('neurons')
title('Wake(J=0.005):stimuli-response')

subplot(2,3,2)
imagesc(Spk_poststim_SO_norm(sig_neuron_ind,:))
hold on
scatter(cenmass_SO(sig_neuron_ind), 1:length(sig_neuron_ind),30,'+w')
xlabel('time from stimuli (ms)'); 
title('SO(J=0.0005):stimuli-response')

subplot(2,3,3)
imagesc(Spk_train_conv_UP_norm(sig_neuron_ind,:))
hold on
scatter(cenmass_UP(sig_neuron_ind), 1:length(sig_neuron_ind),30,'+w')
xlabel('time from UP-onset (ms)'); 
title('SO(J=0.0005):UP-state')

subplot(2,3,4)
scatter(cenmass_W(Ind_E), Ind_E_pos,100,'.r')
hold on
scatter(cenmass_W(Ind_I), Ind_I_pos,100,'.b')
set(gca, 'YDir','reverse');xlim([0 100])

subplot(2,3,5)
scatter(cenmass_SO(Ind_E), Ind_E_pos,100,'.r')
hold on
scatter(cenmass_SO(Ind_I), Ind_I_pos,100,'.b')
set(gca, 'YDir','reverse');xlim([0 100])

subplot(2,3,6)
scatter(cenmass_UP(Ind_E), Ind_E_pos,100,'.r')
hold on
scatter(cenmass_UP(Ind_I), Ind_I_pos,100,'.b')
set(gca, 'YDir','reverse');xlim([0 100])




%% graph theory

clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
%%%%%%%%%%%%%%%%%%%%
load('20201129T161858_Izhi_run.mat') % trend with g
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
addpath('E:\poissonTutorial')

g_peak_store = trend.g_peak_store;

n_permutations   = 100;
N                = param.N;
connect          = initial.connect;
% g_peak           = initial.g_peak ;
N_inh            = round(N*0.2);
N_ex             = N - N_inh;

i_trend = 11
g_peak = squeeze(g_peak_store(i_trend,:,:));

g_peak_strong       = full(g_peak);
g_peak_strong_tmp   = g_peak_strong(find(connect(:,1:N_ex)));
sm_threshold        = prctile(g_peak_strong_tmp,95);
g_peak_strong(g_peak_strong(:,1:N_ex) <= sm_threshold) = 0;
connmat               = (g_peak_strong ~= 0);
connmat(:,(N_ex+1):N) = 0; %%%%%%%%%%%%%%%%%%%%%%%%%% remove I -> E
% connmat((N_ex+1):N,:) = 0; %%%%%%%%%%%%%%%%%%%%%%%%%% remove E -> I
% connmat = (connmat + connmat' > 0); %%% sym

%%%% remove un-used neurons %%%%%%%%
% used_ind = find(sum(connmat + connmat',1) ~= 0);
% connmat = connmat(used_ind, used_ind);

matrix_locs = find(tril(connmat+1-eye(length(connmat))));
nconnections = sum(connmat(:));

% figure
% imagesc(connmat)

clustcoef_real = clustcoef_asym(connmat);
pathlengths_real = pathlen_asym(connmat);

clustercoefficients_random = zeros(1,n_permutations);
pathlengths_random = zeros(1,n_permutations);
parfor permi=1:n_permutations
    disp(['perm ' num2str(permi)])
    
    % generate random network...
    connmat_random = zeros(size(connmat));
    edges2fill     = randsample(matrix_locs,nconnections);
    connmat_random(edges2fill) = 1;
    
    clustercoefficients_random(permi) = clustcoef_asym(connmat_random);
    pathlengths_random(permi) = pathlen_asym(connmat_random);
end

% now compute permuted small-world-network-ness
for permi=1:n_permutations
    whichnetworks2use = randsample(1:n_permutations,2);
    swn_permutations(permi)  = ( clustercoefficients_random(whichnetworks2use(1))/clustercoefficients_random(whichnetworks2use(2)) ) / ( pathlengths_random(whichnetworks2use(1))/pathlengths_random(whichnetworks2use(2)) );
end

swn_real       = ( clustcoef_real/mean(clustercoefficients_random) ) / ( pathlengths_real/mean(pathlengths_random) );
swn_z = (swn_real-mean(swn_permutations(isfinite(swn_permutations))))/std(swn_permutations(isfinite(swn_permutations)));
%
% figure
% histogram(pathlengths_random)
% figure
% histogram(clustercoefficients_random)

figure
[y,x] = hist(swn_permutations,50);
h=bar(x,y,'histc');
set(h,'linestyle','none')
hold on
plot(repmat(swn_real,1,2),get(gca,'ylim')/2,'m','linew',3)
title([ 'Small-world-networkness: Z=' num2str(swn_z) ', p=' num2str(1-normcdf(abs(swn_z))) ])
xlabel('swn value')
ylabel('Count')


[t,s] = find(connmat);
G = digraph(s,t);
figure
plot(G)


%% graph theory, plot graph
clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
addpath('E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\Script\functions_izhi')
addpath('E:\poissonTutorial')
%%%%%%%%%%%%%%%%%%%%
load('20201129T161858_Izhi_run.mat') % trend with g
load('20201201_W_SO_FR_PC_type.mat')

n_permutations   = 100;
N                = param.N;
connect          = initial.connect;
% g_peak           = initial.g_peak ;
N_inh            = round(N*0.2);
N_ex             = N - N_inh;

g_peak_store = trend.g_peak_store;

swn_trend = [];
n_trend = size(g_peak_store,1);

figure
for i_trend = 2:n_trend
    
    %     i_trend = 11
    
    g_peak = squeeze(g_peak_store(i_trend,:,:));
    
    g_peak_strong       = full(g_peak);
    g_peak_strong_tmp   = g_peak_strong(find(connect(:,1:N_ex)));
    sm_threshold        = prctile(g_peak_strong_tmp,95);
    g_peak_strong(g_peak_strong(:,1:N_ex) <= sm_threshold) = 0;
    connmat               = (g_peak_strong ~= 0);
    connmat(:,(N_ex+1):N) = 0; %%%%%%%%%%%%%%%%%%%%%%%%%% remove I -> E
    % connmat((N_ex+1):N,:) = 0; %%%%%%%%%%%%%%%%%%%%%%%%%% remove E -> I
    
    
    %%%% remove un-used neurons %%%%%%%%
    used_ind = find(sum(connmat + connmat',1) ~= 0);
    
    %%%% remove small net from EE EI
    smallnet_EE_EI = [19 398 246 120 421];
    used_ind = setdiff(used_ind,used_ind(smallnet_EE_EI));
    
    connmat = connmat(used_ind, used_ind);
    
    A_N_used   = A_N(used_ind);
    B_N_used   = B_N(used_ind);
    C_N_used   = C_N(used_ind);
    D_N_used   = D_N(used_ind);
    FR_W_used  = FR_W(used_ind);
    FR_SO_used = FR_SO(used_ind);
    PC_W_used  = PC_W(used_ind);
    PC_SO_used = PC_SO(used_ind);
    
    ex_used = (used_ind <= 800);
    in_used = (used_ind > 800);
    
    
    %
    [t,s] = find(connmat);
    G = digraph(s,t);
    
    % figure
    % imagesc(connmat)
    
    % subplot(2,5,i_trend-1)
    % figure
    % p = plot(G,'MarkerSize',5);
    % title(['after ' num2str(i_trend-1) 's training'])
    
    
    [~, clustcoef_real_vec]  = clustcoef_asym(connmat);
    cdata = clustcoef_real_vec;
    
    figure
    subplot(2,2,1)
    p = plot(G,'MarkerSize',3);
    p.NodeCData = FR_W_used ;
    colormap jet;
    h = colorbar;
    ylabel(h, 'FR W')
    subplot(2,2,2)
    p = plot(G,'MarkerSize',3);
    p.NodeCData = FR_SO_used ;
    colormap jet;
    h = colorbar;
    ylabel(h, 'FR SO')
    subplot(2,2,3)
    p = plot(G,'MarkerSize',3);
    p.NodeCData = PC_W_used ;
    colormap jet;
    h = colorbar;
    ylabel(h, 'PC W')
    subplot(2,2,4)
    p = plot(G,'MarkerSize',3);
    p.NodeCData = PC_SO_used ;
    colormap jet;
    h = colorbar;
    ylabel(h, 'PC SO')
    
    
    
    central_type = {'hubs','authorities','outdegree','indegree','betweenness','outcloseness','incloseness','pagerank'};
    for ct = 1:length(central_type)
        
        cdata = centrality(G,central_type{ct});
        
        figure('Renderer', 'painters', 'Position', [10 10 1300 400])
        subplot(2,6,1)
        scatter(cdata(find(ex_used)),FR_W_used(find(ex_used)),'.r')
        ylabel('FR W ex')
        xlabel(central_type{ct})
        lsline
        subplot(2,6,2)
        scatter(cdata(find(ex_used)),FR_SO_used(find(ex_used)),'.r')
        ylabel('FR SO ex')
        xlabel(central_type{ct})
        lsline
        subplot(2,6,7)
        scatter(cdata(find(ex_used)),PC_W_used(find(ex_used)),'.r')
        ylabel('PC W ex')
        xlabel(central_type{ct})
        lsline
        subplot(2,6,8)
        scatter(cdata(find(ex_used)),PC_SO_used(find(ex_used)),'.r')
        ylabel('PC SO ex')
        xlabel(central_type{ct})
        lsline
        
        subplot(2,6,3)
        scatter(cdata(find(in_used)),FR_W_used(find(in_used)),'.b')
        ylabel('FR W in')
        xlabel(central_type{ct})
        lsline
        subplot(2,6,4)
        scatter(cdata(find(in_used)),FR_SO_used(find(in_used)),'.b')
        ylabel('FR SO in')
        xlabel(central_type{ct})
        lsline
        subplot(2,6,9)
        scatter(cdata(find(in_used)),PC_W_used(find(in_used)),'.b')
        ylabel('PC W in')
        xlabel(central_type{ct})
        lsline
        subplot(2,6,10)
        scatter(cdata(find(in_used)),PC_SO_used(find(in_used)),'.b')
        ylabel('PC SO in')
        xlabel(central_type{ct})
        lsline
        
        subplot(2,6,[5 6 11 12])
        p = plot(G,'MarkerSize',3);
        p.NodeCData = cdata ;
        colormap jet;
        h = colorbar;
        ylabel(h, central_type{ct})
        mxx = prctile(cdata,90);
        caxis([min(cdata) mxx]);
        
        
    end
    
    
    
end


%% graph theory, trend swn

clear all
cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
%%%%%%%%%%%%%%%%%%%%
load('20201129T161858_Izhi_run.mat') % trend with g


n_permutations   = 100;
N                = param.N;
connect          = initial.connect;
% g_peak           = initial.g_peak ;
N_inh            = round(N*0.2);
N_ex             = N - N_inh;

g_peak_store      = trend.g_peak_store;

n_trend           = size(g_peak_store,1);
clustcoef_trend   = zeros(n_trend,1);
pathlengths_trend = zeros(n_trend,1);
swn_trend         = zeros(n_trend,1);
for i_trend = 1:n_trend
    
    
    g_peak = squeeze(g_peak_store(i_trend,:,:));
    
    g_peak_strong       = full(g_peak);
    g_peak_strong_tmp   = g_peak_strong(find(connect(:,1:N_ex)));
    sm_threshold        = prctile(g_peak_strong_tmp,95);
    g_peak_strong(g_peak_strong(:,1:N_ex) <= sm_threshold) = 0;
    connmat               = (g_peak_strong ~= 0);
    connmat(:,(N_ex+1):N) = 0; %%%%%%%%%%%%%%%%%%%%%%%%%% remove I -> E
    % connmat((N_ex+1):N,:) = 0; %%%%%%%%%%%%%%%%%%%%%%%%%% remove E -> I
    % connmat = (connmat + connmat' > 0); %%% sym
    
    %%%% remove un-used neurons %%%%%%%%
    % used_ind = find(sum(connmat + connmat',1) ~= 0);
    % connmat = connmat(used_ind, used_ind);
    
    matrix_locs = find(tril(connmat+1-eye(length(connmat))));
    nconnections = sum(connmat(:));
    
    % figure
    % imagesc(connmat)
    
    clustcoef_real = clustcoef_asym(connmat);
    pathlengths_real = pathlen_asym(connmat);
    
    clustercoefficients_random = zeros(1,n_permutations);
    pathlengths_random = zeros(1,n_permutations);
    parfor permi=1:n_permutations
        disp([num2str(i_trend-1) 's run perm ' num2str(permi)])
        
        % generate random network...
        connmat_random = zeros(size(connmat));
        edges2fill     = randsample(matrix_locs,nconnections);
        connmat_random(edges2fill) = 1;
        
        clustercoefficients_random(permi) = clustcoef_asym(connmat_random);
        pathlengths_random(permi) = pathlen_asym(connmat_random);
    end
    
    % now compute permuted small-world-network-ness
    for permi=1:n_permutations
        whichnetworks2use = randsample(1:n_permutations,2);
        swn_permutations(permi)  = ( clustercoefficients_random(whichnetworks2use(1))/clustercoefficients_random(whichnetworks2use(2)) ) / ( pathlengths_random(whichnetworks2use(1))/pathlengths_random(whichnetworks2use(2)) );
    end
    
    swn_real       = ( clustcoef_real/mean(clustercoefficients_random) ) / ( pathlengths_real/mean(pathlengths_random) );
    
    clustcoef_trend(i_trend) = clustcoef_real;
    pathlengths_trend(i_trend) = pathlengths_real;
    swn_trend(i_trend) = swn_real;
end

cd E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run
save([datestr(now, 'yyyymmddTHHMMSS') '_swn_trend' ],'-v7.3');


%%%%%%%
figure
subplot(2,1,1)
yyaxis right
plot(0:10,clustcoef_trend)
ylabel('cluster coefficient')
yyaxis left
plot(0:10,pathlengths_trend)
ylabel('path length')
xlabel('s training')
xlim([0 10])

subplot(2,1,2)
plot(0:10,swn_trend,'k')
xlim([0 10])
xlabel('s training')
ylabel('small-world-netowrkness')
