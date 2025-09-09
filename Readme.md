# What system are we solving?
$$
\begin{cases}
        \dfrac{\partial n_i}{\partial t} + \dfrac{\partial}{\partial z}(n_i u_i) = n_i(\nu_i - \nu_{iw}) \\[2ex] 
        \dfrac{\partial}{\partial t} (n_i u_i) + \dfrac{\partial}{\partial z}(n_i u_i^2) = \dfrac{e}{m_i}n_iE +  n_i(\nu_i u_n - \nu_{iw} u_i) \\[2ex]
        \dfrac{\partial}{\partial t}(\dfrac{3}{2}n_i T_e) + \dfrac{\partial}{\partial z}(\dfrac{5}{2} n_i T_e u_{e} + \kappa_{\perp} \dfrac{\partial T_e}{\partial z}) = -e n_i u_{e} E - S_{coll} - S_{wall} \\[2ex]
        \dfrac{\partial n_n}{\partial t} + u_n \dfrac{\partial n_n}{\partial z} = n_i(\nu_{iw} - \nu_{i})
\end{cases}
$$
where:
$
    E = \dfrac{j}{e n_i \mu_{\perp}} - \dfrac{u_i}{\mu_{\perp}} - \dfrac{1}{n_i}\dfrac{\partial}{\partial z}(n_i T_e)
$

$
    j = \dfrac{U_0 + \int^L_0 (\dfrac{u_i}{\mu_{\perp}} + \dfrac{1}{n_i} \dfrac{\partial }{\partial x} (n_i T_e))dx}{\int_0^L \dfrac{1}{e n_i \mu_{\perp}} dx}
$

$
    \mu_{\perp} = \dfrac{e}{m_e \nu_m} \dfrac{1}{1 + \omega_{ce}^2 / \nu_m^2} 
$

$
    \omega_{ce} = \dfrac{e B }{m_e}
$

$
B = B_0 \exp{\left(-\dfrac{(z - z_0)^2}{2 \delta_B^2}\right)}
$

$
z_0 = 2.5 \; cm, \; B_0 = 150 \; G, \; \delta_{B, in} = 1.1 \; cm, \; \delta_{B, out} = 1.8 \; cm 
$

$
    \nu_m = \nu_{en} + \nu_{ew} + \nu_{ano}
$

$
\nu_{en} = n_n \beta_n 
$

$
\nu_i = n_n \beta_i
$

$\beta_i, \; \beta_n$ - reaction rates from Bolsig+ solver

$
\nu_{iw} = \dfrac{1.2}{R_{out} - R_{in}} \sqrt{\dfrac{T_e}{m_i}} 
$

$
\nu_{ew} = \dfrac{\nu_{iw}}{1 - \delta}  
$

$
\delta = \min(\Gamma(2 + b)a\left(\dfrac{T_e}{e}\right)^b, \; 0.983)
$

$
a = 0.123, \; b = 0.528, \; \Gamma(2 + b) = 1.36
$

$
\nu_{ano} = \alpha_a \omega_{ce}
$

$
\alpha_{a, in} = 1/160, \; \alpha_{a, out} = 1/16
$

$
\kappa_{\perp} = -\dfrac{5}{2} \dfrac{\mu_{\perp}}{e} n_i T_e
$

$
S_{coll} = n_i \nu_i e E_i
$

$
E_i =  12.13 \; eV
$

$
S_{wall} = n_i \nu_{ew} (2 T_e + (1 - \delta) e \varphi_w)
$

$
\varphi_w = \dfrac{T_e}{e}\ln{\left\{ (1 - \delta) \sqrt{\dfrac{m_i}{2 \pi m_e}}\right\}}
$

$
u_n = 10^4 \; cm / s
$



