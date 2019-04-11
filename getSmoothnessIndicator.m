function [smth_idc_arr] = getSmoothnessIndicator(q,w)
  restrict_mat = [ 0.7349768 , 0.37068122,-0.1445425 , 0.03888448;
                  -0.0923266 , 0.5923266 , 0.5923266 ,-0.0923266 ;
                   0.03888448,-0.1445425 , 0.37068122, 0.7349768  ];
  prolong_restrict_mat = [ 1.17382422,-0.23592625, 0.06210202;
                           0.31577941, 0.80735482,-0.12313423;
                          -0.12313423, 0.80735482, 0.31577941;
                           0.06210202,-0.23592625, 1.17382422 ];
  rho_mat = q(:,:,1);
  rho_rest_mat = restrict_mat * rho_mat;
  rho_rest_prol_mat = prolong_restrict_mat * rho_rest_mat;
  rho_diff2_mat = (rho_mat - rho_rest_prol_mat).^2;
  weights_mat = repmat(w,1,size(rho_mat,2));
  rho_diff2_int_arr = sum(weights_mat .* rho_diff2_mat, 1);
  rho_int_arr = sum(weights_mat .* rho_mat, 1);
  smth_idc_arr = log10(rho_diff2_int_arr / rho_int_arr);
end