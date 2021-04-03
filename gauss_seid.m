function x = gauss_seid(A,b,xprev,tol,max_iter)
  n = length(A);
  L_star = zeros(n);
  U = zeros(n);
  
  for i = 1:n
    L_star(i,i) = A(i,i);
    for j = i+1:n
      U(i,j) = A(i,j);
    endfor
  endfor
  
  for i = 2:n
    for j = 1:i-1
      L_star(i,j) = A(i,j);
    endfor
  endfor
  
  epsi = norm(A*xprev-b);
  k = 0; 
  while (epsi > tol) || (k < max_iter)
    x = inv(L_star)*(b-(U*xprev));
    xprev = x;
    epsi = norm(A*xprev-b);
    k = k+1;   
    fprintf('%11.0f %4.10f %4.10f %4.10f %4.10f\n',[k;x] )
  endwhile
endfunction
