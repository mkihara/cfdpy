# Fractional step method

## 計算の流れ

Fractional step method には様々な方法がある。今回は、以下の手順で計算する。

1. 圧力項以外を時間積分する
1. ポアソン方程式を導出し、圧力場を求める
1. 圧力項を時間積分する

## 圧力のポアソン方程式の導出

$$
\frac{\partial \bm{u}}{\partial t} = - \nabla p
$$

上式とクランク・ニコルソン法より下式を得る。

$$
u^{n+1} - u^* = - \frac{\Delta t}{2} \left( \nabla p^{n+1} + \nabla p^n \right)
$$

上式の発散をとる。

$$
\Leftrightarrow \nabla^2 \left( p^{n+1} +p^n \right) = \frac{2}{\Delta t} \left( - \nabla \cdot u^{n+1} + \nabla \cdot u^* \right)
$$

ここで、連続の式より、$\nabla \cdot u^{n+1} = 0$ であり、$p'=\frac{p^{n+1}+p^n}{2}$ とおくと、上式は、

$$
\Leftrightarrow \nabla^2 p' = \frac{1}{\Delta t} \nabla \cdot u^*
$$

となる。上式より、$p^{n+1}$ を得るには、$p^n$ が必要であることがわかる。計算の最初のステップで $p^n$ が不明な場合のためのポアソン方程式も以下で導出しておく。

$$
\frac{\partial \bm{u}}{\partial t} = - \nabla p
$$

と Backward Euler method を用いると、以下が得られる。

$$
u^{n+1} = u^* - \Delta t \nabla p^{n+1}
$$

上式の発散をとると、以下が得られる。

$$
\nabla^2 p^{n+1} = \frac{1}{\Delta t} \left( - \nabla \cdot u^{n+1} + \nabla \cdot u^* \right)
$$

ここで、連続の式より、$\nabla \cdot u^{n+1} = 0$ であるから、上式より、以下のポアソン方程式が得られる。

$$
\nabla^2 p^{n+1} = \frac{1}{\Delta t} \nabla \cdot u^*
$$