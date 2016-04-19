% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT options file
% 
% References
%   - http://www.coin-or.org/Ipopt/documentation/
%   - http://drops.dagstuhl.de/opus/volltexte/2009/2089/pdf/09061.WaechterAndreas.Paper.2089.pdf
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output (C.1)
options.ipopt.print_level                   = 5;
options.ipopt.print_user_options            = 'no';
options.ipopt.print_options_documentation   = 'no';

% Termination (C.2)
options.ipopt.tol                           = 1e-8; % (Opt1)
options.ipopt.dual_inf_tol                  = 1;    % (Opt2)
options.ipopt.constr_viol_tol               = 1e-4; % (Opt3)
options.ipopt.compl_inf_tol                 = 1e-4; % (Opt4), Optimal Solution Found if (Opt1),...,(Opt4) fullfiled

options.ipopt.acceptable_iter               = 3;    % (Acc1)
options.ipopt.acceptable_tol                = 1e10; % (Acc2)
options.ipopt.acceptable_constr_viol_tol    = 1e10; % (Acc3)
options.ipopt.acceptable_dual_inf_tol       = 1e10; % (Acc4)
options.ipopt.acceptable_compl_inf_tol      = 1e10; % (Acc5)
options.ipopt.acceptable_obj_change_tol     = 1e-3; % (Acc6), Solved To Acceptable Level if (Acc1),...,(Acc6) fullfiled

options.ipopt.max_iter                      = 500;
options.ipopt.max_cpu_time                  = 3000;

% Barrier Parameter (C.6)
options.ipopt.mu_strategy = 'adaptive';

% Line Sarch (C.8)
%options.ipopt.accept_every_trial_step = 'yes';
%options.ipopt.line_search_method = 'cg-penalty';

% Restoration Phase (C.10)
%options.ipopt.soft_resto_pderror_reduction_factor = 100;
%options.ipopt.required_infeasibility_reduction    = 0.9999;

% Quasi-Newton (C.13)
options.ipopt.hessian_approximation         = 'limited-memory';
options.ipopt.limited_memory_max_history    = 6;
options.ipopt.limited_memory_initialization = 'scalar2';

% Derivative Test (C.14)
%options.ipopt.derivative_test               = 'first-order';
%options.ipopt.derivative_test_tol           = 1e-3;