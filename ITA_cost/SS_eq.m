function z = SS_eq( alpha1,alpha2,alpha3,z_grid_c,z_grid_d,lambda_grid,lambda_grid_c,lambda_grid_d,im,ih,func_var )
mu_c = func_var.mu_c;
mu_d = func_var.mu_d;
if (im>(alpha1+alpha3)*lambda_grid_c(1)) && (ih > (alpha2+alpha3)*lambda_grid_d(1))
    z=[0,0];
else
    % LHS cash FOC, 1/beta*1/R_z - 1 = i_m 
    % Cash such that liquidity premium is equal to cash interest rate with
    % cash only
    zm_can = interp1((alpha1+alpha3)*lambda_grid_c,(1-mu_c)*z_grid_c,im); 
    % Deposits such that liquidity premium is equal to deposit interest rate
    zh_can = interp1((alpha2+alpha3)*lambda_grid_d,(1-mu_d)*z_grid_d,ih); 
    % Liquidity premium for deposits when holding cash only
    Ah = alpha2*lambda_grid_d(1)+alpha3*interp1((1-mu_c)*z_grid_c,lambda_grid_c,zm_can); 
    % Liquidity premium for cash when holding deposits only
    Am = alpha1*lambda_grid_c(1)+alpha3*interp1((1-mu_d)*z_grid_d,lambda_grid_d,zh_can);    
    if Ah<ih
        z = [zm_can,0];
    elseif Am<im
        z = [0,zh_can];
    else
        % Liquidity premium and amount of cash for type 1 meetings (cash only) when holding
        % both cash and deposits
        AA1 = im-ih+alpha2*lambda_grid_d;
        zm_zh = interp1(alpha1*lambda_grid_c,(1-mu_c)*z_grid_c,AA1);
        ID = ~isnan(zm_zh);
        zm_zh = zm_zh(ID);
        z_grid1 = z_grid_c(ID);
        z_grid_d1 = z_grid_d(ID);
        % Find deposit such that liquidity premium is equal to interest
        % rate
        z_upper = max(max(z_grid_d),max(z_grid_c));
        %AA2 = alpha2*lambda_grid_d(ID)+alpha3*interp1((1-mu_d)*z_grid_d,lambda_grid_d,min(z_grid_d1+zm_zh,z_upper))-ih;
        AA2 = alpha2*lambda_grid_d(ID)+alpha3*interp1((1-mu_d)*z_grid_d+(1-mu_c)*z_grid_c,lambda_grid,min(z_grid_d1+zm_zh,z_upper))-ih;
        [~,ID] = min(abs(AA2));
        zh = z_grid1(ID);
        z = [interp1(z_grid_d1,zm_zh,zh),zh];     
    end 
end
end