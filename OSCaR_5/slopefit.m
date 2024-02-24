function [kp] = slopefit(kp, kpv, F_og, const)

for i=1:length(kpv)
    kp.(kpv{i}).zi=F_og(kp.(kpv{i}).xi,kp.(kpv{i}).yi);
    %kp.(kpv{i}).zf=F_og(kp.(kpv{i}).xf,kp.(kpv{i}).yf);
    for j=1:length(kp.(kpv{i}).zi)
        zfguess1=F_og(kp.(kpv{i}).xf(j),kp.(kpv{i}).yf(j));
        yfguess1=kp.(kpv{i}).yf(j);
        ydel=(kp.(kpv{i}).yf(j)-kp.(kpv{i}).yi(j));
        zdel=(kp.(kpv{i}).zi(j)-zfguess1);
        if abs(zdel)>0.0001
            while abs(zdel)>0.0001
                newspan=-sqrt(abs(ydel^2-(kp.(kpv{i}).zi(j)-zfguess1)^2)); %negative to keep with span convention
                newyf=kp.(kpv{i}).yi(j)+newspan;
                zfguess=F_og(kp.(kpv{i}).xf(j),newyf);
                zdel=zfguess1-zfguess;
                zfguess1=zfguess;
                %ydel=newspan;
            end
            kp.(kpv{i}).zf(j)=zfguess;
            kp.(kpv{i}).yf(j)=newyf;
        else
            kp.(kpv{i}).zf(j)=zfguess1;
            kp.(kpv{i}).yf(j)=kp.(kpv{i}).yi(j)+ydel;
        end
    end
    kp.(kpv{i}).zf=kp.(kpv{i}).zf';
    kp.(kpv{i}).slp=atand((kp.(kpv{i}).zf-kp.(kpv{i}).zi)./kp.(kpv{i}).len);
    kp.(kpv{i}).minlen(kp.(kpv{i}).slp>=const.wp_slope_range | kp.(kpv{i}).slp<=-const.wp_slope_range)=const.wp_slope_outrange;
    kp.(kpv{i}).minlen(kp.(kpv{i}).slp<const.wp_slope_range & kp.(kpv{i}).slp>-const.wp_slope_range)=const.wp_slope_inrange;
    %kp.(kpv{i}).yfa=kp.(kpv{i}).yf+tand(kp.(kpv{i}).slp).*(kp.(kpv{i}).zi-kp.(kpv{i}).zf);

end
end