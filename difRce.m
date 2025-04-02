    % Definice symbolických proměnných
    syms t t0 y0 R1 R2 C U real
    
    % % Definice symbolické funkce y(t)
    syms y(t)
    
    % Rovnice:  -C*y'(t) + y(t)/R2 + (U0 - y(t))/R1 = 0
    ode = -C*diff(y, t) + y/R2 + (U - y)/R1 == 0;
    
    % Počáteční podmínka: y(t0) = y0
    cond = y(t0) == y0;
    
    % Symbolické řešení pomocí dsolve
    sol = dsolve(ode, cond, 'IgnoreAnalyticConstraints', true);
    
    % Zobrazení výsledku v příkazovém okně:
    disp('Symbolické řešení y(t):')
    pretty(sol)