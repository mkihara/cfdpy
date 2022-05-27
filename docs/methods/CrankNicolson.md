# クランク・ニコルソン法 (Crank–Nicolson method)

2次精度の陰解法。

$$
u^{n+1} = u^{n} + \frac{\Delta t}{2} \left( \frac{\partial u^{n+1}}{\partial t} + \frac{\partial u^n}{\partial t} \right) + O(\Delta t^3)
$$

## 導出

$u^n$ と $\frac{\partial u^{n}}{\partial t}$ をそれぞれテイラー展開すると以下が得られる。

$$
u^{n} = u^{n+1} - \Delta t \frac{\partial u^{n+1}}{\partial t} + \frac{1}{2!} \Delta t^2 \frac{\partial^2 u^{n+1}}{\partial t^2} + O(\Delta t^3)
$$

$$
\frac{\partial u^{n}}{\partial t}  = \frac{\partial u^{n+1}}{\partial t} - \Delta t \frac{\partial^2 u^{n+1}}{\partial t^2} + O(\Delta t^2)
$$

2式から $\frac{\partial^2 u^{n+1}}{\partial t^2}$ を消去する。

$$
u^{n} = u^{n+1} - \Delta t \frac{\partial u^{n+1}}{\partial t} + \frac{\Delta t}{2} \left( \frac{\partial u^{n+1}}{\partial t} - \frac{\partial u^n}{\partial t} \right) + O(\Delta t^3)
$$

$$
\Leftrightarrow u^{n} = u^{n+1} - \frac{\Delta t}{2} \left( \frac{\partial u^{n+1}}{\partial t} + \frac{\partial u^n}{\partial t} \right) + O(\Delta t^3)
$$

$$
\Leftrightarrow u^{n} = u^{n+1} - \frac{\Delta t}{2} \left( \frac{\partial u^{n+1}}{\partial t} + \frac{\partial u^n}{\partial t} \right) + O(\Delta t^3)
$$

上式より、以下のクランク・ニコルソン法が得られる。
