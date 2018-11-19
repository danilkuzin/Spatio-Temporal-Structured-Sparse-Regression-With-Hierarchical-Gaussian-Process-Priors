for t = 1 : T
    if (any(any(isnan(Qf(t).var_beta))))
        display(['Qf var beta for ' num2str(t)]);
    end
    if (any(isnan(Qf(t).mean_beta)))
        display(['Qf mean beta for ' num2str(t)]);
    end
    if (isnan(Qf(t).bern_param))
        display(['Qf bern param for ' num2str(t)]);
    end
    
    if (any(any(isnan(Qg(t).var_beta))))
        display(['Qg var beta for ' num2str(t)]);
    end
    if (any(isnan(Qg(t).mean_beta)))
        display(['Qg mean beta for ' num2str(t)]);
    end
    
    if (any(any(isnan(Qh(t).var_gamma))))
        display(['Qh var gamma for ' num2str(t)]);
    end
    if (any(isnan(Qh(t).mean_gamma)))
        display(['Qh mean gamma for ' num2str(t)]);
    end
    if (isnan(Qh(t).bern_param))
        display(['Qh bern param for ' num2str(t)]);
    end
    
    if (any(any(isnan(Qr(t).var_gamma))))
        display(['Qr var gamma for ' num2str(t)]);
    end
    if (any(isnan(Qr(t).mean_gamma)))
        display(['Qr mean gamma for ' num2str(t)]);
    end
    if (any(any(isnan(Qr(t).var_mu))))
        display(['Qr var mu for ' num2str(t)]);
    end
    if (any(isnan(Qr(t).mean_mu)))
        display(['Qr mean mu for ' num2str(t)]);
    end

    if (any(any(isnan(Qu(t).var_mu_left))))
        display(['Qu var mu left for ' num2str(t)]);
    end
    if (any(isnan(Qu(t).mean_mu_left)))
        display(['Qu mean mu left for ' num2str(t)]);
    end
    if (any(any(isnan(Qu(t).var_mu_right))))
        display(['Qu var mu right for ' num2str(t)]);
    end
    if (any(isnan(Qu(t).mean_mu_right)))
        display(['Qu mean mu right for ' num2str(t)]);
    end

end