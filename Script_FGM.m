%        -------------------------------------------------------------
%
%          NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
%
%          Copyright 2005, The Johns Hopkins University
%             School of Medicine and Whiting School of Engineering.
%             All rights reserved.
%             For research use only; commercial use prohibited.
%             Distribution without permission of Joseph L. Greenstein or
%             Raimond L. Winslow not permitted.
%           (jgreenst@jhu.edu, rwinslow@jhu.edu)
%
%          Copyright 2006,  University of California, San Diego
%            Department of Bioengineering
%            All rights reserved.
%            For research use only; commercial use prohibited.
%            Distribution without permission of Sarah N. Flaim or
%            Andrew D. McCulloch not permitted.
%           (shealy@bioeng.ucsd.edu, amcculloch@ucsd.edu)      
%
%          Original model and manuscript:
%          Name of Program: 40-State coupled LCC-RyR Model
%          Version: version 1.0
%          Date: September 2005
%         
%          Greenstein, J. L., R. Hinch, and R.L. Winslow. 2006. Protein Geometry and Placement in the
%          Cardiac Dyad Influence Macroscopic Properties of Calcium-Induced Calcium-Release. 
%          Biophys J. 90(1): 77-91.
%
%          Modified model and manuscript:
%          Description: Canine ventricular myocyte models of endo, epi and
%          M cell ionic current and excitation-contraction coupling.
%          Major modifications: Replacement of Hodgkin-Huxley INa model with Clancy Markov model, including late INa
%                              Inclusion of novel open state inactivation for Ito1 model
%                              Inclusion of selected transmurally varying ion channel parameters                                       
%          Version: version 2.0
%          Date: March 2006
%        
%          Flaim, S.N, Giles, W.R., and McCulloch, A.D 2006. Contributions of sustained INa and IKv4.3 to 
%          transmural heterogeneity of early repolarization and arrhythmogenesis in canine left ventricular 
%          myocytes
%          Am J Physiol Heart Circ Physiol 291: H2617-H2629, 2006
%
%
%        ------------------------------------------------------------- 
 

global Vclamp_flag Period Shift Pulse_duration Pulse_amplitude Nstates_CaRU IKs_array
global Counter T_array ICaL_array JCaL_array JRyR_array CaSSavg_array celltype

Counter = 1;
Length_array = 10000;
options1 = odeset('MaxStep',0.1,'Stats','on');       %   Set integrator options

T_array = zeros(Length_array,1);        %   Initialize large arrays
T_array(1) = 123;                       %   Large value here prevents Counter from incrementing on first call to fcn_40state
ICaL_array = zeros(Length_array,1);
JCaL_array = zeros(Length_array,1);
JRyR_array = zeros(Length_array,1);
CaSSavg_array = zeros(Length_array,1);
IKs_array = zeros(Length_array,1);
T = [];
Y = [];

Vclamp_flag = 0;                        %   = 0 for APs, = 1 for voltage clamp 
Period = 2000;                          %   Period between stimuli or voltage clamp (ms)
Shift=10;                               %   Delay between start of Period and stimulus (ms)
Pulse_duration=5;                       %   Duration of stimulus (ms)
Pulse_amplitude=-10.0;                  %   Amplitude of stimulus current (A/F)
Num_beats = 15;                         %   Simulation runs on interval [0, Num_beats*Period]

Nstates_CaRU = 40;                      %   Number of states in the coupled LCC-RyR model
celltype = 'endo'                        %   Cell subtype = epi/mid/endo 
                                        
if Vclamp_flag                          %   Set initial conditions of states
    IC_CaRU = zeros(1,40);              %   LCC/RyR states
    IC_CaRU(1) = 1;
    IC_global_VC = [ -100.0100    0.0001    0.9995    0.9995   10.0000  134.5668    0.0001    0.7500    0.1033    0.9781 ...
                        0.0001    0.9683    0.0134    0.0001    0.0000    0.0000    0.0153    0.0027    0.0002    0.0000 ...
                        0.0000    0.8232    0.0522    0.0012    0.0000    0.0000    0.1192    0.0034    0.0007    0.0001 ...
                        0.0000    0.9998    0.0002    0.0000    0.0000    0.0000    0.9995]; %   Ionic model states
    IC_VC = [ IC_CaRU IC_global_VC 0];
    IC = IC_VC;
    IC(78) = 0.0;  %ICs for INa Markov model
    IC(42) = 5.40077e-10; %UO
    IC(43) = 9.4146e-09; %UIF
    IC(44) = 0.00183429; %UC2
    IC(77) = 0.990169; %UC3 
    IC(79) = 2.39346e-12; %UIM1
    IC(81) = 0.00765092; %IC3
    IC(82) = 1.41733e-05; %IC2
    IC(83) = 3.20257e-17; %UIM2
    IC(84) = 0.000330056; %LC3
    IC(85) = 6.11428e-07; %LC2
    IC(86) = 4.0614e-10; %LC1
    IC(87) = 1.80026e-13; %LO
    IC(80) = 1 - sum([IC(2:4),IC(77),IC(79:80)]);
else 
    IC = load('NewerICs/IC_endo_2000.txt');
end

if Vclamp_flag              %   Uses Vclamp_flag to choose simulation(s) below
    AP_flag = 0;
    Gain_flag = 1;
else
    AP_flag = 1;
    Gain_flag = 0;
end    
                       
if AP_flag                  % Run APs.  Requires Vclamp_flag = 0
    for i = 1:Num_beats
        T_range = [(i-1)*Period (i-1)*Period+Shift]
        [Ta,Ya] = ode23t(@fcn_fgm,T_range,IC,options1);
        T_range = [(i-1)*Period+Shift i*Period]
        IC = Ya(end,:);
        [Tb,Yb] = ode23t(@fcn_fgm,T_range,IC);
        T = [T; Ta; Tb];
        Y = [Y; Ya; Yb];
        IC = Yb(end,:);
    end
    figure(1); hold on; plot(T,Y(:,Nstates_CaRU+1),'k-')
    IC = Y(end,:);
end

if Gain_flag                % Voltage clamp
    VC = [ -20.01:10:-10.01 0.01:10:70.01];
    clear JCaLp JRyRp
    rVC = round(VC);
    amplitudes = zeros(length(VC),1);
    IKs = [];
    T = [];
    for i = 1:length(VC)       
        VC(i)
        T_range = [0 100];
        IC(Nstates_CaRU+1) = -50.01;
        [T1,Y1] = ode23t(@fcn_fgm,T_range,IC);
        T_range = [100 3100];
        IC = Y1(end,:);
        IC(Nstates_CaRU+1) = 50.01;
        [T2,Y2] = ode23t(@fcn_fgm,T_range,IC);
        T_range = [3100 3400];
        IC = Y2(end,:);
        IC(Nstates_CaRU+1) = -25.01;;
        [T3,Y3] = ode23t(@fcn_fgm,T_range,IC,options1);
        T = [T;T_array];
        T_array = T_array(1:Counter);
        IKs_array = IKs_array(1:Counter,:);
        JCaL_array = JCaL_array(1:Counter);
        ICaL_array = ICaL_array(1:Counter);
        figure(1);hold on;plot(T_array,IKs_array(:,1),'r-');
        figure(2);hold on;plot(T_array,ICaL_array(:,1),'r-');
    end
end

T_array = T_array(1:Counter);           %   Truncate unused space from these storage arrays
ICaL_array = ICaL_array(1:Counter);
JCaL_array = JCaL_array(1:Counter);
JRyR_array = JRyR_array(1:Counter);
CaSSavg_array = CaSSavg_array(1:Counter);
IKs_array = IKs_array(1:Counter,:);



