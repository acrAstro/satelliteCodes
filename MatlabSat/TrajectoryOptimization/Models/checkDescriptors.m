function [CD,DDI,DD2] = checkDescriptors(CD,DDI,DD2)

if strcmpi(CD,'Classical') == 1
    CD = 'Classical';
elseif strcmpi(CD,'Nonsingular') == 1
    CD = 'Classical';
else
    error('Please check Chief Descriptor\n');
end

if strcmpi(DDI,'Cartesian') == 1
    DDI = 'Cartesian';
elseif strcmpi(DDI,'Relative Classical')
    DDI = 'Relative Classical';
elseif strcmpi(DDI,'Relative Nonsingular')
    DDI = 'Relative Nonsingular';
else
    error('Please check Deputy Descriptor');
end

if strcmpi(DD2,'Cartesian') == 1
    DD2 = 'Cartesian';
elseif strcmpi(DD2,'Relative Classical')
    DD2 = 'Relative Classical';
elseif strcmpi(DD2,'Relative Nonsingular')
    DD2 = 'Relative Nonsingular';
else
    error('Please check Deputy Descriptor');
end

end