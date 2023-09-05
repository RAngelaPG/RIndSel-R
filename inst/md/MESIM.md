## Molecular Eigen Selection Index Method (MESIM).

MESIM is a generalization of ESIM to the case in which is incorporated MM information to IS, similarly as in the IS Lande and Thompson (1990). Following the basic idea of Kempthorne and Nordskog (1959), Ceron-Rojas {et al}. (2008b), maximized the correlation between $Y_M ={\rm{\bf{\beta} '}}_p{\rm{\bf p}} +{\rm{\bf{\beta}'}}_s{\rm{\bf s}} ={\rm{\bf{\beta}'}}_M{\rm{\bf p}}_{ps}$  and $Z_M ={\rm{\bf{\theta } '}}_1{\rm{\bf g}} +{\rm{\bf{\theta}'}}_2{\rm{\bf s}} ={\rm{\bf {\theta}'}}_M{\rm{\bf g}}_{gs}$  , $\rho_{Y_M Z_M}^2$ , where ${\rm{\bf p}$ and ${\rm{\bf g}}$ have been defined in Equations 1 and 2, and ${\rm {\bf s}}$ is the vector of records of the additive effects of QTLs 
associated with MM, ${\rm{\bf{p} '}}_{ps}=\left [{{\begin{array}{*{20} c} {{\rm{\bf{p} '}}} \hfill &{{\rm{\bf{s}'}}} \hfill \\ \end{array}}} \right]$, ${\rm{\bf{g} '}}_{gs}=\left [ {{\begin{array}{*{20} c} {{\rm{\bf{g} '}}} \hfill &{{\rm{\bf{s}'}}} \hfill \\ \end{array}}} \right]$, ${\rm{\bf{\beta} '}}_M = \left [ {{\begin{array}{*{20} c} {{\rm{\bf{\beta} '}}_P} \hfill &{{\rm{\bf{\beta}'}}_s} \hfill \\
\end{array}}} \right]$, and ${\rm{\bf{\theta} '}}_M = \left [ {{\begin{array}{*{20} c} 
{{\rm{\bf{\theta} '}}_1} \hfill &{{\rm{\bf{\theta}'}}_2} \hfill \\ \end{array}}} \right]$. In the illustrative example of Lande and Thompson (1990), $Y_M = \beta_p*p + \beta_s*s$, and $Z_M = \theta_p*G_P + \theta_s*s$ , where $p$ denotes the plant height and $s = x_1 \alpha_1 + x_2 \alpha_2 + ... + x_5 \alpha_5$. 

Again, because $\rho_{Y_M Z_M }^2$  is invariant to changes in is maximized scale $\rho_{Y_M Z_M ^2}$ under the constraints ${\rm{\bf{\beta} '}}_M{\rm{\bf T \beta}}_M = 1$ $\beta '_M*T*\beta_M = 1$ and ${\rm{\bf{\theta} '}}_M{\rm{\bf K \theta}}_M =1$ $\theta '_M*K*\theta_M = 1$, so in MESIM, it is necessary maximize 
\[
\Phi = \left ({{\rm{\bf{\theta} '}}_M{\rm{\bf K \beta}}_M} \right) ^ 2 - \mu 
\left ({{\rm{\bf{\beta} '}}_M{\rm{\bf T \beta}}_M -1} \right) - \omega 
\left ({{\rm{\bf{\theta} '}}_M{\rm{\bf K \theta}_M} -1} \right) 
\] 
With respect to ${\rm{\bf \beta}}_M $, ${\rm{\bf \theta}}_M$, $\mu$ , and $\omega$, where ${\rm{\bf \beta}}_M$ is the vector of coefficients $\mbox{MESIM}$, ${\rm{\bf \theta}}_M$ is the vector of coefficients $Z_M ={\rm{\bf{\theta} '}}_M{\rm{\bf g}}_{gs}$ and $\mu$ and $\omega$ are Lagrange multipliers. Ceron-Rojas {et al}. (2008b) found that the solution is equal 
\[ 
\label{EQ4} 
({\rm{\bf Q}} - \mu{\rm{\bf I}}){\rm{\bf \beta}}_M ={\rm{\bf 0}}, 
\]
where ${\rm{\bf Q}}={\rm{\bf T}} ^{-1}{\rm{\bf K}}$. Thus, in $\mbox{MESIM}$, the value that maximizes $\rho_{Y_M Z_M} ^ 2$ under the 
constraints ${\rm{\bf{\beta} '}}_M{\rm{\bf T \beta}}_M = 1$ and ${\rm {\bf{\theta} '}}_M{\rm{\bf K \theta}}_M = 1$ is the first eigenvalue ($\mu$) Matrix ${\rm{\bf Q}}$, and the vector for building $Y_M$ (With maximum correlation with $Z_M ={\rm{\bf{\theta}'}}_M{\rm{\bf g}}_{gs})$ is the first eigenvector (${\rm{\bf \beta}}_M$) matrix ${\rm{\bf Q}}$.
