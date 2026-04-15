using Random
using Printf

Random.seed!(1234)

# Glosario de variables:
# eta - viscosidad
# rho - densidad
# sigma - posición 
# E - campo eléctrico
# q - carga
# r - radio

#constantes
const rho_aceite = 1064.0
const g = 9.80665
const e = 1.602176634e-19
const d = 6.0e-3 # distancia entre placas
const rho_aire_base = 1.225
const eta_base = 1.86e-5

#funciones para la simulación
function dvdt(v,r,rho_aire,eta,m,sigma,q,E)
    F_g = (4/3) * π * r^3 * (rho_aceite-rho_aire) * g
    F_d=6*pi*eta*r*v
    F_e=q*E*sigma
    dvdt = (F_e + F_g - F_d) / m
    return dvdt
end

# Método de Runge-Kutta de 4to orden para resolver la ecuación diferencial
function runge_kutta(v0,r,rho_aire,eta,m,sigma,q,E,dt)
    k1 = dvdt(v0,r,rho_aire,eta,m,sigma,q,E)
    k2 = dvdt(v0 + 0.5 * dt * k1,r,rho_aire,eta,m,sigma,q,E)
    k3 = dvdt(v0 + 0.5 * dt * k2,r,rho_aire,eta,m,sigma,q,E)
    k4 = dvdt(v0 + dt * k3,r,rho_aire,eta,m,sigma,q,E)
    v_next = v0 + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
    return v_next
end

# Simulación de la caída y subida de la gota
function millikan(r,rho_aire,eta,m,sigma,q,E,dt,t_max)
    n_steps = Int(ceil(t_max / dt))
    v = zeros(n_steps)
    v[1] = 0.0
    
    for i in 1:(n_steps-1)
        v[i+1] = runge_kutta(v[i], r, rho_aire, eta, m, sigma, q, E, dt)
    end
    return v
end

N = 30
resultados = []

println("Gota    r (μm)   n   vf (m/s)   vs (m/s)   q_est (C)   q_real (C)   err %")
println("-" ^ 95)

# Simulación para N gotas
for i in 1:N
    r = rand() * 4e-6 + 1e-6 
    n_elec = rand(1:8)  # Carga aleatoria (entre 1 y 8 electrones)
    q = n_elec * e

    m = (4/3) * π * r^3 * rho_aceite

    # Variaciones para simular condiciones experimentales
    eta = eta_base * (1 + 0.01*(rand()-0.5)) 
    rho_aire = rho_aire_base * (1 + 0.01*(rand()-0.5))
    sigma = 1.0 + 0.02*(rand()-0.5)
    E = 5e4 * (1 + 0.02*(rand()-0.5))

    # Simulación sin campo eléctrico
    v_caida = millikan(r, rho_aire, eta, m, sigma, q, 0.0, 1e-6, 5.0)
    v_caida_final = v_caida[end] 

    # Simulación con campo eléctrico
    v = millikan(r, rho_aire, eta, m, sigma, q, E, 1e-6, 5.0)
    v_final = v[end]

    # Estimación del radio y carga a partir de la velocidad de caída
    r_est = sqrt((9 * eta * v_caida_final) / (2 * g * (rho_aceite - rho_aire)))
    F_g = (4/3) * π * r_est^3 * (rho_aceite - rho_aire) * g

    m_est = (4/3) * π * r_est^3 * rho_aceite

    # Simulación de la subida con el campo eléctrico
    v_subida = millikan(r_est, rho_aire, eta, m_est, sigma, q, -E, 1e-6, 5.0)
    v_subida_final = v_subida[end]

    # Estimación de la carga a partir de la velocidad de subida
    q_est = (6 * π * eta * r_est * (v_caida_final - v_subida_final)) / E
    q_real = q

    # Cálculo del error porcentual
    q_est_entre_e = q_est / e # Estimación del número de electrones a partir de q_est
    n_elec_float = Float64(n_elec) # Convertir n_elec a Float64 para el cálculo del error porcentual
    err_percent = abs(q_est_entre_e - n_elec_float) / n_elec_float * 100 # Error porcentual basado en el número de electrones estimado vs real
    
    # Conversión del radio a micrómetros para la presentación de resultados
    r_um = r * 1e6

    # Impresión de resultados
    @printf("%2d    %.3f   %d   %.6e   %.6e   %.3e   %.3e   %.2f %%\n", 
        i, r_um, n_elec, v_final, v_subida_final, q_est, q_real, err_percent)

    push!(resultados, (
    gota=i, 
    r=r_um, 
    n=n_elec, 
    vf=v_final, 
    vs=v_subida_final, 
    q_est=q_est, 
    q_real=q_real, 
    err=err_percent
    ))

end