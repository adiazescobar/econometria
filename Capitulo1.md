---
date: 2023-07-21
title: "Notas de Clase Econometría"
---

# Ejercicios de Montecarlo {#ejercicios-de-montecarlo .unnumbered}

La tarea de la econometría involucra la estimación de parámetros, como por ejemplo, la media poblacional o los coeficientes de una regresión lineal, a partir de una muestra de datos extraída del mundo real. No obstante, no sólo nos interesa la estimación en sí, sino también la precisión de esta, es decir, cuán próxima se encuentra nuestra predicción de los valores reales. Cabe mencionar que un estimador es, en su esencia, una función de variables aleatorias y, por ende, se comporta como una variable aleatoria.

En la realidad, nos encontramos con la posibilidad de observar una única muestra de un tamaño predefinido n, que nos permite estimar un solo parámetro. Sin embargo, los experimentos de Monte Carlo nos brindan la oportunidad de replicar esta realidad cuantas veces deseemos, digamos R repeticiones. En cada una de estas réplicas, obtenemos una muestra distinta de tamaño n de la población original. Esto nos da la capacidad de calcular los parámetros de interés en múltiples ocasiones, obteniendo en cada intento un resultado diferente. La distribución empírica de estas estimaciones sucesivas tiende a aproximarse al valor verdadero del estimador.

Para facilitar la comprensión, recurramos a un ejemplo proporcionado por Kennedy (2006). Supongamos que disponemos de 45 observaciones de las variables $y$ y $x$. Tenemos conocimiento de que $y$ se ha generado a partir de las observaciones de $x$, siguiendo la fórmula: $y=\beta x+\epsilon$, donde $\epsilon$ representa el error con una media de cero y una varianza de $\sigma^{2}$.


Es importante destacar que esta relación lineal no incluye una constante o término independiente. Por consiguiente, sólo contamos con un parámetro desconocido, $\beta$, que es la pendiente de $x$. Esto implica que $y_{1}$, la primera observación de $y$, se generó multiplicando $x_{1}$, la primera observación de $x$, por el parámetro desconocido $\beta$ y sumando luego $\epsilon_{1}$, el primer error. Este último se obtiene seleccionando aleatoriamente un error de un conjunto de errores con media cero y varianza $\sigma^{2}$. El procedimiento para generar las 44 observaciones restantes es idéntico.

Nuestro interés radica en estimar el parámetro desconocido $\beta$. Dos sugerencias se presentan: el profesor A propone la fórmula $\hat\beta_A=\sum y/\sum x$ mientras que el profesor B sugiere $\hat\beta_B=\sum xy/\sum x^{2}$. Supongamos que, con los datos que poseemos, calculamos los siguientes coeficientes: $\hat\beta_A=2.34$ y $\hat\beta_B=2.73$. ¿Cuál de estos dos parámetros sería su elección?

Volviendo a nuestro proceso generador de datos $y = \beta x + \epsilon$
y reemplazándolo en cada una de las sugerencias obtenemos lo siguiente:

$$\hat\beta_A=\frac{\sum(\beta x+\epsilon)}{\sum x} = \beta +\frac{\sum\epsilon}{\sum x}$$

$$\hat\beta_B=\frac{\sum(\beta x+\epsilon)}{\sum x^{2}} = \beta +\frac{\sum x \epsilon}{\sum x^{2}}$$

De estas formulas podemos observar que los parámetros estimados son equivalentes
al $\beta$ verdadero incrementado por una expresión que incluye los valores
desconocidos de $\epsilon$. Recordemos que $\epsilon$ tiene una media
igual a cero, lo que nos lleva a esperar que sea negativo aproximadamente
la mitad de las veces y positivo la otra mitad. En otras palabras, el
segundo término en cada una de las fórmulas tiende a ser un valor cercano a cero, lo que sugiere que las fórmulas propuestas por los
profesores A y B son razonables: ambas ofrecen una estimación adecuada
de $\beta$.

Consideremos la estimación propuesta por el profesor A: $\hat\beta_A=\sum(y)/\sum(x)=2.34$. La cercanía entre $\hat\beta_A$ y $\beta$ depende intrínsecamente de los errores seleccionados aleatoriamente del conjunto de errores para generar $y$. Si predominan errores positivos y grandes al recopilar los datos, entonces $\hat\beta_A$ sobrestimaría $\beta$. Por el contrario, si predominan los errores negativos, $\hat\beta_A$ subestimaría $\beta$. Por tanto, al seleccionar un conjunto típico de errores, es probable que $\hat\beta_A$ se acerque mucho a $\beta$. No obstante, es importante destacar que el conjunto de 45 errores seleccionados determina el estimador producido por la fórmula de $\hat\beta_A$, por lo que resulta imposible determinar qué tan cerca se encuentra la estimación $\hat\beta_A=2.34$ del valor verdadero.

Supongamos por un momento que conocemos el valor verdadero de $\beta$.

$$Y= \beta x+\epsilon$$

Adoptemos que $\beta=2$ y que $\epsilon$ sigue una distribución normal con media cero y varianza cuatro, o sea, $\epsilon\sim N(0,4)$.

Podemos sustituir $\beta$ en las fórmulas propuestas por los profesores A y B para obtener:

$$\hat{\beta}_{A}=2+\frac{\sum\epsilon}{\sum x}$$

$$\hat{\beta}_{B}=2+\frac{\sum x\epsilon}{\sum x^{2}}$$

En Stata, podemos calcular $\hat{\beta}_{A}$ y $\hat{\beta}_{B}$ de la siguiente manera:



```stata
set obs 45 
set seed 1111 
generate x = runiform()*10 
generate epsilon = rnormal(0,2)
sum epsilon 
Variable | Obs Mean Std. Dev. Min Max 
epsilon | 45 -.1525543 2.101546 -3.9534 3.599413 
generate y = 2*x+epsilon
gen unos = 1
mata 
y = st_data(., "y") 
x = st_data(., "x") 
unos = st_data(., "unos")
beta_a = (unos'x)^(-1)*unos'y 
beta_a  //output: 1.966063212 
beta_b = (x'x)^(-1)*x'y 
beta_b  //output: 1.94818922 
end
```

Es evidente que, en este caso, la mayoría de las perturbaciones seleccionadas aleatoriamente fueron negativas. Por lo tanto, tanto $\hat{\beta}_{A}$ como $\hat{\beta}_{B}$ tienden a subestimar a $\beta$.

Los resultados son intrínsecamente dependientes de las 45 observaciones de $\epsilon$ seleccionadas de manera aleatoria. Si repitiéramos el experimento, es probable que los valores de ambos estimadores variaran. Para simplificar este proceso, podríamos implementar un programa en Stata que nos permita repetir el experimento cuantas veces deseemos, seleccionando cada vez una muestra aleatoria de perturbaciones.


```stata
capture program drop betas 
program define betas, rclass 
capture drop epsilon y 
generate epsilon = rnormal(0,2) 
generate y = 2*x+epsilon 
mata: 
y = st_data(., "y") 
x = st_data(., "x") 
unos = st_data(., "unos") 
epsilon = st_data(., "epsilon") 
beta_a = (unos'x)^(-1)*unos'y 
st_numscalar("beta_a" , beta_a) 
beta_b = (x'x)^(-1)*x'y 
st_numscalar("beta_b" , beta_b) 
end
return scalar beta_a = beta_a 
return scalar beta_b = beta_b 
di "Beta_A = " beta_a 
di "Beta_B = " beta_b 
end
```
Si lo llamamos una vez, obtenemos los siguientes resultados:

```stata
betas 
Beta_A = 2.0175583 Beta_B = 2.0523701
````
Al parecer las perturbaciones escogidas aleatoriamente se encontraban
balanceadas entre valores negativos y positivos. Intentemos una vez más:

```stata
. betas Beta_A = 1.9541299 Beta_B = 1.9513717
```

Una vez más las perturbaciones escogidas aleatoriamente fueron en su
mayoría negativas. Una vez más:

```stata
. betas Beta_A = 1.9591604 Beta_B = 2.0594701
```

Y así sucesivamente, digamos 2000 veces. En lugar de generar 2000
repeticiones de $\epsilon$ podemos utilizar el comando simulate de Stata
de la siguiente manera:

```stata
simulate beta_hat_a=r(beta_a) beta_hat_b=r(beta_b), reps(2000) nodots:
betas command: betas beta_hat_a: r(beta_a) beta_hat_b: r(beta_b)
```

Este programa nos produce 2000 valores hipotéticos de $\hat{\beta}_{A}$
y de $\hat{\beta}_{B}$. Podemos realizar un summarize para ver los
estadísticos de este estimador:

```stata
sum beta_hat_A beta_hat_B

Variable | Obs Mean Std. Dev. Min Max
beta_hat_A | 2000 1.999181 .0668024 1.699073 2.301055 
beta_hat_B | 2000 1.999335 .0530487 1.787391 2.238478
```

Tal vez es más fácil realizar el análisis con un gráfico:

::: stlog
twoway (kdensity beta_hat_A, lwidth(medthick)) /// \> (kdensity
beta_hat_B, lwidth(medthick)), ytitle(Densidad) xtitle(Betas) /// \>
xline(2, lwidth(medthick) lpattern(dash) lcolor(black))
title(Distribución \> Muestral 2000 Repeticiones) /// \> legend(order(1
\"Beta A\" 2 \"Beta B\") box fcolor(none) lcolor(none) cols(1))
:::

::: center
[Distribución Muestral Estimadores A y B]{.image .placeholder
original-image-src="stata/cap1" original-image-title="fig:"}
:::

Como se puede ver claramente en el gráfico los valores hipotéticos
$\hat{\beta}_{A}$ y de $\hat{\beta}_{B}$, nunca toman valores muy
lejanos de $\beta$, ya que necesitaríamos muestras muy atípicas de las
perturbaciones.

Estos histogramas estiman la distribución de valores posibles de
$\hat\beta_A$ y $\hat\beta_B$. Esta distribución se conoce como la
distribución muestral de cada uno de los parámetros. Recuerde que un
estimador es insesgado si su distribución muestral se encuentra centrada
alrededor del verdadero parámetro a ser estimado. En este caso, la
distribución muestral de los dos estimadores cumple este criterio y por
lo tanto podríamos decir que los dos estimadores son insesgados. Ahora
bien, entre estas dos distribuciones, **¿cuál considera usted que es la
mejor?**

Veamos otro ejemplo, también tomado de Kennedy (2006), en el cual
incluimos una constante: $y=\alpha+\beta x+\epsilon$. Si usamos la
definición del Profesor B, $\hat\beta=\sum(xy)/\sum(x^2)$ e incluimos
$y$ definido anteriormente, tenemos que:

$$\hat\beta_B=\beta +\sum x \alpha / \sum x^2 +\sum x \epsilon / \sum x^2$$

Como se mencionó anteriormente el último término de esta ecuación debe
ser cero, por lo tanto en este caso $\hat\beta$, al contrario que en el
ejemplo anterior, está sesgado por una cantidad igual a
$\sum x \alpha / \sum x^2$.

El siguiente estimador, visto por ustedes en clase de econometría, es
sugerido por el profesor C:

$$\hat\beta_c=\sum (x-\bar{x})(y-\bar{y})/\sum (x-\bar{x})^2$$
**Demuestre que el estimador $\hat\beta_C$ es insesgado**

Si nosotros fuéramos a escoger el estimador basándonos en el criterio
del estimador menos sesgado, deberíamos escoger $\hat\beta_C$ en lugar
de $\hat\beta_B$. Pero antes de saber cuál de los dos estimadores es
mejor realicemos el mismo ejercicio de Monte Carlo para obtener las
distribuciones muestrales y así decidir cuál de los dos estimadores
tiene mejores propiedades.

En Stata podemos modificar el programa anterior donde asumimos que el
proceso generador de datos es:

$$\begin{aligned}
y & = & \alpha+\beta x+\epsilon\\
& = & 1+2x+\epsilon\quad\textrm{donde }\epsilon\sim N(0,4)
\end{aligned}$$

Es decir asumimos que $\alpha=1$, que $\beta=2$ y que
$\epsilon\sim N(0,4)$

Podemos hacer un simple programa que nos permita realizar múltiples
repeticiones:

::: stlog
. capture program drop betas . program define betas, rclass 1. capture
drop epsilon y 2. generate epsilon = rnormal(0,2) 3. generate y =
1+2\*x+epsilon 4. mata: y = st_data(., \"y\") 5. mata: x = st_data(.,
\"x\") 6. mata: unos = st_data(., \"unos\") 7. mata: x2 = st_data(.,
\"unos x\") 8. mata: epsilon= st_data(., \"epsilon\") 9. mata: beta_a =
(unos'x)(-1)\*unos'y 10. mata: st_numscalar(\"beta_a\" , beta_a) 11.
mata: beta_b = (x'x)-1\*x'y 12. mata: st_numscalar(\"beta_b\" , beta_b)
13. mata: beta_c =invsym(x2'x2)\*x2'y 14. mata: st_numscalar(\"beta_c\"
, beta_c\[2,1\]) 15. return scalar beta_a = beta_a 16. return scalar
beta_b = beta_b 17. return scalar beta_c = beta_c 18. di \"Beta_A = \"
beta_a 19. di \"Beta_B = \" beta_b 20. di \"Beta_C = \" beta_c 21. end .
. betas Beta_A = 2.1794083 Beta_B = 2.1456709 Beta_C = 2.0761487
:::

Al igual que en el caso anterior, este programa nos calcula los
estimadores propuestos por cada uno de los profesores. Usando el comando
**simulate** podemos obtener 2000 valores hipotéticos de
$\hat{\beta}_{B}$ y de $\hat{\beta}_{C}$, los cuales pueden emplearse
para construir la distribución muestral de los dos estimadores:

::: stlog
. sum beta_hat_B beta_hat_C Variable Obs Mean Std. Dev. Min Max
beta_hat_B 2000 2.149365 .05416 1.974667 2.336887 beta_hat_C 2000
2.000123 .0960512 1.686663 2.320208
:::

Gráficamente tenemos

::: center
[Distribución Muestral Estimadores B y C]{.image .placeholder
original-image-src="stata/cap1b" original-image-title="fig:"}
:::

En la gráfica se puede ver claramente que el estimador $\hat{\beta}_{B}$
es un estimador sesgado, mientras que el estimador $\hat{\beta}_{C}$ es
insesgado. Sin embargo, $\hat{\beta}_{B}$ tiene menor varianza que
$\hat{\beta}_{C}$. Dadas estas características ¿cuál considera usted que
es el mejor estimador de $\beta$?

Recordemos que buscamos estimadores que sean insesgados y que sean
eficientes. Un estimador es eficiente si su varianza es la menor entre
todas las varianzas de los estimadores insesgados. Al reconocer la
eficiencia como una propiedad deseable en un estimador estamos dejando
excluidos a estimadores sesgados incluso cuando su varianza sea muy
pequeña, como es el caso de $\hat{\beta}_{B}$. En ocasiones podríamos
tolerar cierto sesgo a cambio de una varianza reducida. Un criterio que
nos permite tener en cuenta las situaciones el las que tenemos
estimadores con sesgo tolerable unido a una varianza pequeña es el
criterio del Error Cuadrático Medio, el cual se define como:

$$\begin{aligned}
ECM(\beta) & = & E(\hat{\beta}-\beta)^{2}\\
 & = & E(\hat{\beta}^{2})-2\beta E(\hat{\beta})+\beta^{2}\\
 & = & [E(\hat{\beta}^{2})-E(\hat{\beta})^{2}]+[E(\hat{\beta})^{2}-2\beta E(\hat{\beta})+\beta^{2}]\\
 & = & [E(\hat{\beta}^{2})-E(\hat{\beta})^{2}]+[E(\hat{\beta})-\beta)]^{2}\\
 & = & Var(\hat{\beta})+Sesgo(\hat{\beta})^{2}
\end{aligned}$$

En este caso elegimos el estimador que minimice el ECM. En otras
palabras, elegimos el estimador que cometa, en promedio, el menor error
en la estimación. Es decir que tenga menos sesgo y menor varianza. Por
ejemplo, si estamos comparando dos estimadores insesgados, según este
criterio escogeríamos aquellos estimadores que tengan varianza mínima.
Si comparamos un estimador insesgado y otro sesgado no necesariamente
debemos escoger el insesgado.

Siguiendo con nuestro ejemplo podemos calcular el ECM para cada
estimador de la siguiente manera:

::: stlog
. gen ECM_B = (beta_hat_B - 2)2 . gen ECM_C = (beta_hat_C - 2)2 . sum
ECM_B ECM_C Variable Obs Mean Std. Dev. Min Max ECM_B 2000 .0252416
.0167352 8.36e-06 .1134926 ECM_C 2000 .0092212 .0132326 1.61e-09
.1025331
:::

En este caso el estimador con menor ECM es$\hat{\beta}_{C}$. Usando la
última formula de la ecuación anterior podemos descomponer el error
cuadrático medio como la suma de dos componentes: la varianza del
estimador más su sesgo al cuadrado.

::: stlog
. di \"Varianza Beta_B= \" varbb Varianza Beta_B= .00293184 . di \"Sesgo
Beta_B 2= \" sesgo_B2 Sesgo Beta_B 2= .02230977 . di \"ECM_B= \"
varbb+sesgo_B2 ECM_B= .02524161 . . di \"Varianza Beta_C= \" varbc
Varianza Beta_C= .00922121 . di \"Sesgo Beta_C 2= \" sesgo_C2 Sesgo
Beta_C 2= 1.503e-08 . di \"ECM_C= \" varbc+sesgo_C2 ECM_C= .00922123
:::

# Ejercicios

1.  Demuestre que el estimador $\hat{\beta}_{C}$ es insesgado.

2.  Basado en el do file, de la clase de Stata, realice lo siguiente:

    1.  Genere 100 valores de $x$ de una distribución uniforme entre 2 y
        8.

    2.  Genere 100 valores de $\epsilon$ de una distribución normal
        estándar. Multiplique cada uno de estos valores por 3 y
        súmele 2. ¿Cómo se distribuye $\epsilon$?

    3.  Cree la variable $y$ usando la siguiente fórmula:
        $y=\alpha+\beta x+\epsilon$, donde $\alpha=2.5$, $\beta=5$ y
        $\epsilon$ tiene las propiedades enunciadas en el punto
        anterior.

    4.  Calcule $\hat{\beta}_{B}$ y $\hat{\beta}_{C}$ usando las
        definiciones del Profesor B y C (ver notas de clase)

    5.  Repita el proceso descrito anteriormente comenzando en el
        numeral (b) hasta que generé 5000 valores de $\hat{\beta}_{B}$ y
        $\hat{\beta}_{C}$ (pista: use el comando simulate y cambie el
        número de repeticiones). Calcule la media y la varianza e
        interprete.

    6.  Haga un gráfico de la distribución muestral de ambos
        estimadores. Interprete.

    7.  Calcule el Error Cuadrático Medio. Interprete.
