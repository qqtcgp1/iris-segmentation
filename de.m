function [ optimum, best_value ] = de( S_struct )
%#codegen
%********************************************************************
% Script file for the initialization and run of the differential
% evolution optimizer.
%********************************************************************

if ~isfield(S_struct, 'F_VTR')
    % F_VTR		"Value To Reach" (stop when ofunc < F_VTR)
    S_struct.F_VTR = 0.00001;
end

if ~isfield(S_struct, 'I_D')
    % I_D		number of parameters of the objective function
    S_struct.I_D = 2;
end

% FVr_minbound,FVr_maxbound   vector of lower and bounds of initial population
%    		the algorithm seems to work especially well if [FVr_minbound,FVr_maxbound]
%    		covers the region where the global minimum is expected
%               *** note: these are no bound constraints!! ***
if ~isfield(S_struct, 'FVr_minbound')
    S_struct.FVr_minbound = -100*ones(1,S_struct.I_D);
    S_struct.FVr_maxbound = +100*ones(1,S_struct.I_D);
end

if ~isfield(S_struct, 'I_bnd_constr')
    S_struct.I_bnd_constr = 0;  %1: use bounds as bound constraints, 0: no bound constraints
end

if ~isfield(S_struct, 'I_NP')
    % I_NP            number of population members
    S_struct.I_NP = max(10*S_struct.I_D,20);
end

if ~isfield(S_struct, 'I_itermax')
    % I_itermax       maximum number of iterations (generations)
    S_struct.I_itermax = 50000;
end

if ~isfield(S_struct, 'F_weight')
    % F_weight        DE-stepsize F_weight ex [0, 2]
    S_struct.F_weight = 0.80;
end

if ~isfield(S_struct, 'F_CR')
    % F_CR            crossover probabililty constant ex [0, 1]
    S_struct.F_CR = 1;
end

if ~isfield(S_struct, 'I_strategy')
    % I_strategy     1 --> DE/rand/1:
    %                      the classical version of DE.
    %                2 --> DE/local-to-best/1:
    %                      a version which has been used by quite a number
    %                      of scientists. Attempts a balance between robustness
    %                      and fast convergence.
    %                3 --> DE/best/1 with jitter:
    %                      taylored for small population sizes and fast convergence.
    %                      Dimensionality should not be too high.
    %                4 --> DE/rand/1 with per-vector-dither:
    %                      Classical DE with dither to become even more robust.
    %                5 --> DE/rand/1 with per-generation-dither:
    %                      Classical DE with dither to become even more robust.
    %                      Choosing F_weight = 0.3 is a good start here.
    %                6 --> DE/rand/1 either-or-algorithm:
    %                      Alternates between differential mutation and three-point-
    %                      recombination.
    
    S_struct.I_strategy = 1;
end

if ~isfield(S_struct, 'I_refresh')
    % I_refresh     intermediate output will be produced after "I_refresh"
    %               iterations. No intermediate output will be produced
    %               if I_refresh is < 1
    S_struct.I_refresh = 0;
end

if ~isfield(S_struct, 'I_plotting')
    % I_plotting    Will use plotting if set to 1. Will skip plotting otherwise.
    S_struct.I_plotting = 0;
end

%***************************************************************************
% Problem dependent but constant values. For speed reasons these values are
% defined here. Otherwise we have to redefine them again and again in the
% cost function or pass a large amount of parameters values.
%***************************************************************************

if ~isfield(S_struct, 'FVr_bound')
    %-----bound at x-value +/- 1.2--------------------------------------------
    S_struct.FVr_bound = repmat([-100, 100], [1, S_struct.I_D]);
end

%-----Definition of tolerance scheme--------------------------------------
%-----The scheme is sampled at I_lentol points----------------------------
if ~isfield(S_struct, 'I_lentol')
    S_struct.I_lentol   = min(max(S_struct.I_NP*0.5,30), S_struct.I_NP);
end

if ~isfield(S_struct, 'FVr_x')
    S_struct.FVr_x      = linspace(-1,1,S_struct.I_lentol); %ordinate running from -1 to +1
end
if ~isfield(S_struct, 'FVr_lim_up')
    S_struct.FVr_lim_up = ones(1,S_struct.I_lentol);   %upper limit is 1
end
if ~isfield(S_struct, 'FVr_lim_lo')
    S_struct.FVr_lim_lo = -ones(1,S_struct.I_lentol);  %lower limit is -1
end
%-----tie all important values to a structure that can be passed along----
%-----do not change anything below----------------------------------------
% S_struct.FVr_bound  = FVr_bound;
% S_struct.I_lentol   = I_lentol;
% S_struct.FVr_x      = FVr_x;
% S_struct.FVr_lim_up = FVr_lim_up;
% S_struct.FVr_lim_lo = FVr_lim_lo;
%
% S_struct.I_NP         = I_NP;
% S_struct.F_weight     = F_weight;
% S_struct.F_CR         = F_CR;
% S_struct.I_D          = I_D;
% S_struct.FVr_minbound = FVr_minbound;
% S_struct.FVr_maxbound = FVr_maxbound;
% S_struct.I_bnd_constr = I_bnd_constr;
% S_struct.I_itermax    = I_itermax;
% S_struct.F_VTR        = F_VTR;
% S_struct.I_strategy   = I_strategy;
% S_struct.I_refresh    = I_refresh;
% S_struct.I_plotting   = I_plotting;

%-----values just needed for plotting-------------------------------------
% if (I_plotting == 1)
%    FVr_xplot                = linspace(-1.3,1.3,100); %just needed for plotting
%    FVr_lim_lo_plot          = FVr_bound(I_D)-(FVr_bound(I_D)+1)*step(FVr_xplot+1.2)+(FVr_bound(I_D)+1)*step(FVr_xplot-1.2);
%    S_struct.FVr_xplot       = FVr_xplot;
%    S_struct.FVr_lim_lo_plot = FVr_lim_lo_plot;
% end
%********************************************************************
% Start of optimization
%********************************************************************
f = @(x,y) (struct('I_nc', 0, 'FVr_ca', 0, 'I_no', 1, 'FVr_oa', S_struct.objective(x)));
[FVr_x,S_y,~] = deopt(f, S_struct);
optimum = FVr_x;
best_value = S_y.FVr_oa;


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:         S_MSE= objfun(FVr_temp, S_struct)
% Author:           Rainer Storn
% Description:      Implements the cost function to be minimized.
% Parameters:       FVr_temp     (I)    Paramter vector
%                   S_Struct     (I)    Contains a variety of parameters.
%                                       For details see Rundeopt.m
% Return value:     S_MSE.I_nc   (O)    Number of constraints
%                   S_MSE.FVr_ca (O)    Constraint values. 0 means the constraints
%                                       are met. Values > 0 measure the distance
%                                       to a particular constraint.
%                   S_MSE.I_no   (O)    Number of objectives.
%                   S_MSE.FVr_oa (O)    Objective function values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
