function intpoint_time = evaluate_time(cell,r)

tic;
intersectionpoint(cell,r,1);
intersectionpoint(cell,r,2);

intpoint_time = toc;