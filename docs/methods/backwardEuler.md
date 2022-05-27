# 後退オイラー法 (Backward Euler method) 

1次精度の陰解法。

$$
u^{n+1} = u^{n} + \Delta t \frac{\partial u^{n+1}}{\partial t} + O(\Delta t^2)
$$

## 導出

$u^n$ をテイラー展開することで得られる。

$$
u^{n} = u^{n+1} - \Delta t \frac{\partial u^{n+1}}{\partial t} + O(\Delta t^2)
$$
