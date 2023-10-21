    function x=mu_law(y,Xmax,mu)
        x=((exp(log(1+mu)*abs(y)/Xmax)-1)*Xmax/mu)*sign(y);
        x=x*sign(y);
    end