function [CD,DD] = checkDescriptors(CD,DD)

if strcmpi(CD,'Classical') == 1
    CD = 'Classical';
elseif strcmpi(CD,'Nonsingular') == 1
    CD = 'Classical';
else
    error('Please check Chief Descriptor\n');
end

if strcmpi(DD,'Cartesian') == 1
    DD = 'Cartesian';
elseif strcmpi(DD,'Relative Classical')
    DD = 'Relative Classical';
elseif strcmpi(DD,'Relative Nonsingular')
    DD = 'Relative Nonsingular';
else
    error('Please check Deputy Descriptor');
end

end