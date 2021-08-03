% Only teh clinical tatget volume of the segmented patient prostaete is
% taken here.
% Following a frequently prescribed planing dose of 15 Gy (Med Phys 39.5, pp.
% 2904–2929), the objectives are updated in the following way.


cst{6,6}{1}.parameters =  {15};
cst{6,6}{1}.className  =  'DoseObjectives.matRad_SquaredUnderdosing';
cst{8,6}{1}.parameters =  {7.5};
cst{9,6}{1}.parameters =  {7.5};
cst{1,6}{1}.parameters =  {7.5};

cst{5,6}{1}            =  cst{6,6}{1};
cst{5,3}               =  'TARGET';
cst{6,6}               =  [];
cst{6,3}               =  'OAR';

% the planing target was changed to the clinical segmentation of the
% prostate bed.
% In this example, the lymph nodes will not be part of the treatment:
cst{7,6}               =  [];
cst{7,3}               =  'OAR';