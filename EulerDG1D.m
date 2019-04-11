function [q] = EulerDG1D(x,q,h,m,N,CFL,gamma,FinalTime)
% function [q] = EulerDG1D(x,q,h,m,N,CFL,gamma,FinalTime)
% Purpose  : Integrate 1D Euler equation until FinalTime using a DG
%            scheme and 3rd order SSP-RK method
% Initialize operators at Legendre Gauss Lobatto grid
[r,w] = LegendreGL(m); V  = VandermondeDG(m, r); D = DmatrixDG(m, r, V);
Ma = inv(V*V'); S = Ma*D; iV = inv(V);

limiter_type = 2; % MINMOD-TVB: 1; WENO: 2;
shock_detector = 4; % MINMAX: 4;
eps_zero = 1E-12;
ids=[];

% Compute operator for WENO smoothness evaluator
[qW,Xm,Xp] = WENODGWeights(m,iV);

% Initialize extraction vector
VtoE = zeros(2,N);
for j=1:N
  VtoE(1,j) = (j-1)*(m+1)+1; VtoE(2,j) = j*(m+1);
end

% Compute smallest spatial scale timestep
rLGLmin = abs(r(1)-r(2)); 
time = 0; tstep = 0; 

tstep_monitor = 1;

% integrate scheme
while (time<FinalTime)
  if(mod(tstep,tstep_monitor) == 0)
    disp(tstep); disp(time);
    % plot(x,q(:,:,1),'o');
  end
  % Set timestep
  p = (gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1)); 
  c = sqrt(gamma*p./q(:,:,1));
  c = real(c);
  maxvel = max(max(c+abs(q(:,:,2)./q(:,:,1)))); k = CFL*h*rLGLmin/maxvel;
  if (time+k>FinalTime) k = FinalTime-time; end
  
  % Stage 1 of SSPRK
  rhsq  = EulerDGrhs1D(x,q,h,k,m,N,gamma,S,Ma,VtoE,maxvel); 
  q1 = q + k*rhsq;
  
  % Limit solution through characteristics
  [qc,R] = EulerQtoRDG(q1,gamma,V,iV);
  if(limiter_type==1)
    R(:,:,1) = SlopeLimitCSDG(x,R(:,:,1),m,h,N,V,iV);
    R(:,:,2) = SlopeLimitCSDG(x,R(:,:,2),m,h,N,V,iV);
    R(:,:,3) = SlopeLimitCSDG(x,R(:,:,3),m,h,N,V,iV);
  elseif(limiter_type==2)
    if(shock_detector==4)
      phi_rho_arr = getShockDetector4(R(:,:,1),w,h,VtoE);
      phi_rhoE_arr = getShockDetector4(R(:,:,3),w,h,VtoE);
      phi_rho_arr = min(phi_rho_arr,[],1);
      phi_rhoE_arr = min(phi_rhoE_arr,[],1);
      ids = (1.0-phi_rho_arr>eps_zero) .* (1.0-phi_rhoE_arr>eps_zero);
    end
    R(:,:,1) = WENOlimitDG(x,R(:,:,1),m,h,N,V,iV,qW,Xm,Xp,shock_detector,ids);
    R(:,:,2) = WENOlimitDG(x,R(:,:,2),m,h,N,V,iV,qW,Xm,Xp,shock_detector,ids);
    R(:,:,3) = WENOlimitDG(x,R(:,:,3),m,h,N,V,iV,qW,Xm,Xp,shock_detector,ids);
  end
  q1 = EulerRtoQDG(R,qc,gamma,V,iV);
  
  % Stage 2 of SSPRK
  rhsq  = EulerDGrhs1D(x,q1,h,k,m,N,gamma,S,Ma,VtoE,maxvel); 
  q2 = (3*q + q1 + k*rhsq)/4;

  % Limit solution through characteristics
  [qc,R] = EulerQtoRDG(q2,gamma,V,iV);
  if(limiter_type==1)
    R(:,:,1) = SlopeLimitCSDG(x,R(:,:,1),m,h,N,V,iV);
    R(:,:,2) = SlopeLimitCSDG(x,R(:,:,2),m,h,N,V,iV);
    R(:,:,3) = SlopeLimitCSDG(x,R(:,:,3),m,h,N,V,iV);
  elseif(limiter_type==2)
    if(shock_detector==4)
      phi_rho_arr = getShockDetector4(R(:,:,1),w,h,VtoE);
      phi_rhoE_arr = getShockDetector4(R(:,:,3),w,h,VtoE);
      phi_rho_arr = min(phi_rho_arr,[],1);
      phi_rhoE_arr = min(phi_rhoE_arr,[],1);
      ids = (1.0-phi_rho_arr>eps_zero) .* (1.0-phi_rhoE_arr>eps_zero);
    end
    R(:,:,1) = WENOlimitDG(x,R(:,:,1),m,h,N,V,iV,qW,Xm,Xp,shock_detector,ids);
    R(:,:,2) = WENOlimitDG(x,R(:,:,2),m,h,N,V,iV,qW,Xm,Xp,shock_detector,ids);
    R(:,:,3) = WENOlimitDG(x,R(:,:,3),m,h,N,V,iV,qW,Xm,Xp,shock_detector,ids);
  end
  q2 = EulerRtoQDG(R,qc,gamma,V,iV);
  
  % Stage 3 of SSPRK
  rhsq  = EulerDGrhs1D(x,q2,h,k,m,N,gamma,S,Ma,VtoE,maxvel); 
  q = (q + 2*q2 + 2*k*rhsq)/3;

  % Limit solution through characteristics
  [qc,R] = EulerQtoRDG(q,gamma,V,iV);
  if(limiter_type==1)
    R(:,:,1) = SlopeLimitCSDG(x,R(:,:,1),m,h,N,V,iV);
    R(:,:,2) = SlopeLimitCSDG(x,R(:,:,2),m,h,N,V,iV);
    R(:,:,3) = SlopeLimitCSDG(x,R(:,:,3),m,h,N,V,iV);
  elseif(limiter_type==2)
    if(shock_detector==4)
      phi_rho_arr = getShockDetector4(R(:,:,1),w,h,VtoE);
      phi_rhoE_arr = getShockDetector4(R(:,:,3),w,h,VtoE);
      phi_rho_arr = min(phi_rho_arr,[],1);
      phi_rhoE_arr = min(phi_rhoE_arr,[],1);
      ids = (1.0-phi_rho_arr>eps_zero) .* (1.0-phi_rhoE_arr>eps_zero);
    end
    R(:,:,1) = WENOlimitDG(x,R(:,:,1),m,h,N,V,iV,qW,Xm,Xp,shock_detector,ids);
    R(:,:,2) = WENOlimitDG(x,R(:,:,2),m,h,N,V,iV,qW,Xm,Xp,shock_detector,ids);
    R(:,:,3) = WENOlimitDG(x,R(:,:,3),m,h,N,V,iV,qW,Xm,Xp,shock_detector,ids);
  end
  q = EulerRtoQDG(R,qc,gamma,V,iV);
 
  time = time+k; tstep = tstep+1; 
end
return