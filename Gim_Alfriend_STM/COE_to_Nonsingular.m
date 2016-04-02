function nsElems = COE_to_Nonsingular(kepElems,tol)
[~, f] = keplerSolve(kepElems(2), kepElems(6), tol);
nsElems(1) = kepElems(1);
nsElems(2) = kepElems(5) + f;
nsElems(3) = kepElems(3);
nsElems(4) = kepElems(2)*cos(kepElems(5));
nsElems(5) = kepElems(2)*sin(kepElems(5));
nsElems(6) = kepElems(4);
end