function [delta_v_1, delta_v_2, delta_t, a_trasf, e_trasf] = ...
    bitangentTransfer (a_i, e_i, a_f, e_f, type, mu)

%     Trasferimento bitangente per orbite ellittiche
%     ---------------------------------------------------------------------
%     Argomenti di input
%
%     a_i         [km]          asse semi-maggiore iniziale
%     e_i         [-]           eccentricità iniziale
%     a_f         [km]          asse semi-maggiore finale
%     e_f         [-]           eccentricità finale
%     type        [char]        tipo di manovra:
%                               type = 'pa'    pericentro a apocentro
%                               type = 'ap'    apocentro a pericentro
%                               type = 'pp'    pericentro a pericentro
%                               type = 'aa'    apocentro a apocentro
%     mu          [km^3 s^-2]   parametro gravitazionale 
%     ---------------------------------------------------------------------
%     Argomenti di output:
%
%     delta_v_1   [km s^-1]     primo impulso della manovra
%     delta_v_2   [km s^-1]     secondo impulso della manovra
%     delta_t     [s]           tempo di trasferimento
%     ---------------------------------------------------------------------

rp_i = a_i * (1 - e_i);
ra_i = a_i * (1 + e_i);
rp_f = a_f * (1 - e_f);
ra_f = a_f * (1 + e_f);

switch type
   case 'pa'
       % Dati dell'orbita di trasferimento
       rp_trasf = rp_i;
       ra_trasf = ra_f;
       a_trasf = (rp_trasf + ra_trasf)/2;
       e_trasf = abs((ra_trasf - rp_trasf)/(ra_trasf + rp_trasf));
       % Delta-v richiesto
       vp_trasf = sqrt( mu * ( 2/rp_trasf - 1/(a_trasf)));
       vp_i = sqrt( mu * ( 2/rp_i - 1/(a_i)));
       va_trasf = sqrt( mu * ( 2/ra_trasf - 1/(a_trasf)));
       va_f = sqrt( mu * ( 2/ra_f - 1/(a_f)));
       delta_v_1 = abs(vp_trasf - vp_i);
       delta_v_2 = abs(va_f - va_trasf);
       % Tempo di trasferimento richiesto
       delta_t = pi * sqrt(a_trasf^3/mu);
   case 'ap'
       % Dati dell'orbita di trasferimento
       rp_trasf = rp_f;
       ra_trasf = ra_i;
       a_trasf = (rp_trasf + ra_trasf)/2;
       e_trasf = abs((ra_trasf - rp_trasf)/(ra_trasf + rp_trasf));
       % Delta-v richiesto
       va_trasf = sqrt( mu * ( 2/ra_trasf - 1/(a_trasf)));
       va_i = sqrt( mu * ( 2/ra_i - 1/(a_i)));
       vp_trasf = sqrt( mu * ( 2/rp_trasf - 1/(a_trasf)));
       vp_f = sqrt( mu * ( 2/rp_f - 1/(a_f)));
       delta_v_1 = abs(va_trasf - va_i);
       delta_v_2 = abs(vp_f - vp_trasf);
       % Tempo di trasferimento richiesto
       delta_t = pi*sqrt(a_trasf^3 / mu);
   case 'pp'
       % Dati dell'orbita di trasferimento
       rp_trasf = rp_i;
       ra_trasf = rp_f;
       a_trasf = (rp_trasf + ra_trasf)/2;
       e_trasf = abs((ra_trasf - rp_trasf)/(ra_trasf + rp_trasf));
       % Delta-v richiesto
       vp_trasf = sqrt( mu * ( 2/rp_trasf - 1/(a_trasf)));
       vp_i = sqrt( mu * ( 2/rp_i - 1/(a_i)));
       va_trasf = sqrt( mu * ( 2/ra_trasf - 1/(a_trasf)));
       vp_f = sqrt( mu * ( 2/rp_f - 1/(a_f)));
       delta_v_1 = abs(vp_trasf - vp_i);
       delta_v_2 = abs(vp_f - va_trasf);
       % Tempo di trasferimento richiesto
       delta_t = pi*sqrt(a_trasf^3 / mu);  
   otherwise 'aa';
       % Dati dell'orbita di trasferimento
       rp_trasf = ra_i;
       ra_trasf = ra_f;
       a_trasf = (rp_trasf + ra_trasf)/2;
       e_trasf = abs((ra_trasf - rp_trasf)/(ra_trasf + rp_trasf));
       % Delta-v richiesto
       vp_trasf = sqrt( mu * ( 2/rp_trasf - 1/(a_trasf)));
       va_i = sqrt( mu * ( 2/ra_i - 1/(a_i)));
       va_trasf = sqrt( mu * ( 2/ra_trasf - 1/(a_trasf)));
       va_f = sqrt( mu * ( 2/ra_f - 1/(a_f)));
       delta_v_1 = abs(vp_trasf - va_i);
       delta_v_2 = abs(va_f - va_trasf);
       % Tempo di trasferimento richiesto
       delta_t = pi*sqrt(a_trasf^3 / mu);
end


