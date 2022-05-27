# 低記憶容量ルンゲ・クッタ法 (Low-storage Runge-Kutta methods)

$N$変数の数値積分において、一般の$n$次精度ルンゲ・クッタ法のメモリ使用量が$nN$であるのに対して、この手法{cite}`WILLIAMSON198048`。は$2N$で済むことから2N-methodとも呼ばれる。

$j = 1, n$

$$
q_j = a_j q_{j-1} + h f(x_{j-1}) \\
x_j = x_{j-1} + b_j q_j
$$

ここで、$x, f$ は常微分方程式 $\dot x = f(x)$ の変数と関数、$a_j, b_j$は定数である。

```{bibliography}
```
