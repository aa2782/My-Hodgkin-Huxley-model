clear

C = 1; 
E_l = -52 ;
E_Na = 55;
E_k = -75;
g_L = 0.3;
g_Na = 120;
g_k = 36;

dt = 0.01;
t = 0:dt:5000;

V = zeros(1,length(t));
m = zeros(1,length(t));
n = zeros(1,length(t));
h = zeros(1,length(t));

hea = zeros(1,length(t));
hea(floor(2000/dt):floor(4000/dt))=1;

Iapp = 1; 


m_inf = @(V) 1/(1+exp(-(V+40)/9));
h_inf = @(V) 1/(1+exp((V+62)/10));
n_inf = @(V) 1/(1+exp(-(V+53)/16));

tau_m =  0.3;
tau_h = @(V) 1+11/(1+ exp((V+62)/10)); 
tau_n = @(V) 1+6/(1+ exp((V+53)/16));

V(1) = -70; 
m(1) = m_inf(-60);
n(1) = n_inf(-60);
h(1) = h_inf(-60);


for i = 1:length(t)-1
    km1 = (m_inf(V(i))-m(i))/tau_m;
    mav1 = m(i) + dt*km1;
    km2 = (m_inf(V(i))-mav1)/tau_m;
    m(i+1) = m(i)+(km1+km2)*dt/2;

    kn1 = (n_inf(V(i))-n(i))/tau_n(V(i));
    nav1 = n(i) + dt*kn1;
    kn2 = (n_inf(V(i))-nav1)/tau_n(V(i));
    n(i+1) = n(i)+(kn1+kn2)*dt/2;

    kh1 = (h_inf(V(i))-h(i))/tau_h(V(i));
    hav1 = h(i)+dt*kh1;
    kh2 = (h_inf(V(i))-hav1)/tau_h(V(i));
    h(i+1) = h(i)+(kh1+kh2)*dt/2;
    
    k1 = (Iapp*hea(i)-g_k*(n(i)^4)*(V(i)-E_k)-g_Na*(m(i)^3)*h(i)*(V(i)-E_Na)-g_L*(V(i)-E_l))/C;
    av1 = V(i)+dt*k1;
    k2 = (Iapp*hea(i+1)-g_k*(n(i)^4)*(av1-E_k)-g_Na*(m(i)^3)*h(i)*(av1-E_Na)-g_L*(av1-E_l))/C;
    V(i+1) = V(i)+(k1+k2)*dt/2;  
end

plot(t,V, 'LineWidth',1.5);


