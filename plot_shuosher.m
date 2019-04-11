import_weno;
ref_x_arr=opweno5SSP404dt1(:,2);
ref_rho_arr=opweno5SSP404dt1(:,3);
hold on; plot(ref_x_arr,ref_rho_arr,'-','LineWidth',2);

q=real(q);
plot(x,q(:,:,1),'o');

out_data=zeros(size(q,1)*size(q,2),size(q,3)+1);
out_data(:,1)=x(:);
q_2 = reshape(q,[size(q,1)*size(q,2),size(q,3)]);
out_data(:,2)=q_2(:,1);
out_data(:,3)=q_2(:,2);
out_data(:,4)=q_2(:,3);
% dlmwrite('matlab_ShuOsher_P3_NC400_RK3_CFL0.1_MINMODTVB-M100.dat',out_data,'delimiter',',');
% dlmwrite('matlab_ShuOsher_P3_NC400_RK3_CFL0.1_WENO-M10.dat',out_data,'delimiter',',');
