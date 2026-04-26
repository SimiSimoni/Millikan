using Printf

# Constantes
const E_0 = 8.854e-12   # Permisividad del vacío (F/m)
const E_r = 1.0         # Aire ~1

# Parámetros
d = 6.0e-3             # distancia entre placas (m)
A = 5.0e-4             # área de las placas (m²)
C = E_0 * E_r * A / d    # capacitancia
R = 1.0e6              # resistencia (Ω)
V_fuente = 500.0       # voltaje de la fuente (V)
T = R * C              # constante de tiempo

println("Constante de tiempo T = ", T, " s")

# Simulación temporal
dt = 1e-3              # paso de tiempo (s)
t_max = 0.5            # duración de la simulación (s)
n_steps = Int(t_max/dt)

#inicialización de arrays para almacenar resultados
t = collect(0:dt:t_max) # vector de tiempo
Vc = zeros(length(t))  # voltajes 
E = zeros(length(t))   # campos eléctricos
Vdist = []             # potencial en diferentes tiempos

for i in 1:length(t)
    # Carga del capacitor (es la formula que encontre, si encontraron otra se puede cambiar solo esta linea sin problema)
    Vc[i] = V_fuente * (1 - exp(-t[i]/T))
    E[i] = Vc[i] / d

    # Guardar la distribución de potencial en tiempos específicos
    if t[i] in [0.01, 0.05, 0.1, 0.2, 0.5]
        x = range(0, d, length=50)  # posición entre placas
        Vx = Vc[i] .* (x ./ d)      # distribución lineal de potencial
        push!(Vdist, (tiempo=t[i], x=x, Vx=Vx))
    end
end