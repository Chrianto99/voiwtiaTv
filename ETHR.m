data = table2array(oikismoi);
d1 = data(:,1);
d2 = data(:,2);
LOS1 = data(:,3);
LOS2 = data(:,4);
s = data(:,5);
h = data(:,6);
repeater = data(:,7);
service = data(:,8);

cable_loss = [8 * 0.498 + 6 * 0.498 + 1 ,5 * 0.57 + 6 * 0.498 + 1 ,0 * 0.498 + 6 * 0.498 + 1];
EIRP = [38,42,35];
f = 0.51*10^9;
u = 299792458;
wl = 100*u/f;
fr = 10.2*10^9;
wlr = 100*u/fr;

path_loss = NaN*zeros(length(data),1);
repeater_loss = zeros(length(data),1);
Jv = zeros(length(data),1);
Margin = 3*ones(length(data),1);
PE = 5*ones(length(data),1);
EIRP_l = NaN*zeros(length(data),1);
EIRP2 = 22*ones(length(data),1);
receiver_gain = 17*ones(length(data),1);
receiver_gain_hf = 39.7*ones(length(data),1);
cable_loss_l = NaN*zeros(length(data),1);
receiver_Power = NaN*zeros(length(data),1);

for i = 1:length(data)
    
    if (LOS1(i) == 0 ) && (LOS2(i) == 0)
        if (repeater(i) == 0)
            path_loss(i) = 122 + 20*log10(d1(i)) - 20*log10(wl);
            v = h(i)*sqrt(2/wl * (1/s(i) + 1/(d1(i)-s(i))));
            Jv(i) = 6.9 + 20*log10(sqrt((v - 0.1)^2 + 1) + v - 0.1) ;
        else
            path_loss(i) = 122 + 20*log10(d1(i)) - 20*log10(wlr);
            repeater_loss(i) = 122 + 20*log10(d2(i)) - 20*log10(wl);
        end

    elseif (LOS1(i) == 1 ) && (LOS2(i) == 1)
        if (~isnan(h(i)))
            v = h(i)*sqrt(2/wl * (1/s(i) + 1/(d1(i)-s(i))));
            Jv(i) = 6.9 + 20*log10(sqrt((v - 0.1)^2 + 1) + v - 0.1) ;
        end
        path_loss(i) = min(122 + 20*log10(d1(i)) - 20*log10(wl),122 + 20*log10(d2(i)) - 20*log10(wl));
        path_loss(i) = 122 + 20*log10(d1(i)) - 20*log10(wl);
    else 
         path_loss(i) = 122 + 20*log10(d1(i)) - 20*log10(wl);
        if (~isnan(h(i)))
            v = h(i)*sqrt(2/wl * (1/s(i) + 1/(d1(i)-s(i))));
            Jv(i) = 6.9 + 20*log10(sqrt((v - 0.1)^2 + 1) + v - 0.1) ;
        end
        
    end
        EIRP_l(i) = EIRP(service(i));
        cable_loss_l(i) = cable_loss(service(i));
        receiver_Power(i) = EIRP(service(i))  - (path_loss(i) + Jv(i) + PE(i) + Margin(i)  ) + 30;
        receiver_Power_cable(i) = EIRP(service(i))  - (path_loss(i) + Jv(i) + PE(i) + Margin(i) + cable_loss(service(i))) ;
    
end


Total_loss = path_loss + Jv + Margin;
T = table(path_loss,Jv,repeater_loss,Margin,Total_loss,'RowNames',oikismoi1)

P = table(receiver_Power,'RowNames',oikismoi1)

final_Power = receiver_Power_cable' + receiver_gain + 30;
FP = table(EIRP_l,path_loss,Jv,Margin,PE,receiver_gain,cable_loss_l,final_Power,'RowNames',oikismoi1)

final_Power_r1 = EIRP_l - path_loss - Jv - Margin  - 1 + receiver_gain_hf + 30;
RP1 = table(EIRP_l,path_loss,Jv,Margin,receiver_gain_hf,final_Power_r1,'RowNames',oikismoi1)


final_Power_r2 = EIRP2 - repeater_loss - Jv - Margin - PE - 9  + receiver_gain + 30;
RP2 = table(EIRP2,repeater_loss,Jv,Margin,PE,receiver_gain,final_Power_r2,'RowNames',oikismoi1)


