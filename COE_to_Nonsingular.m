function NSElems = COE_to_Nonsingular(eleosc,tol)
[~, f] = keplerSolve(eleosc(2), eleosc(6), tol);
NSElems(1) = eleosc(1);
NSElems(2) = eleosc(5)+f;
NSElems(3) = eleosc(3);
NSElems(4) = eleosc(2)*cos(eleosc(5));
NSElems(5) = eleosc(2)*sin(eleosc(5));
NSElems(6) = eleosc(4);
end